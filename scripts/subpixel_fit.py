"""Fit an optical root tip to original grayscale pixels and measure its chord."""

import numpy as np
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import gaussian_filter1d, map_coordinates
from scipy.optimize import least_squares
from scipy.special import erf


DEFAULT_LENGTH_MM = 0.200


def _intensity(gray, points):
    return map_coordinates(
        gray, [points[:, 1], points[:, 0]], order=3, mode="nearest", prefilter=True
    )


def _blurred_bar(parameters, offsets):
    background, slope, amplitude, center, half_width, blur = parameters
    scale = np.sqrt(2) * blur
    return background + slope * offsets + amplitude * 0.5 * (
        erf((offsets - center + half_width) / scale)
        - erf((offsets - center - half_width) / scale)
    )


def _fit_profile(gray, point, tangent):
    normal = np.array([-tangent[1], tangent[0]])
    offsets = np.arange(-16, 16.001, 0.5)
    values = _intensity(gray, point + offsets[:, None] * normal)
    edges = np.r_[values[:8], values[-8:]]
    background = np.median(edges)
    amplitude = max(3.0, np.percentile(values[np.abs(offsets) <= 8], 90) - background)
    noise = max(0.7, np.median(np.abs(np.diff(edges))))
    excess = np.maximum(values - background, 0) * (np.abs(offsets) < 9)
    center = np.sum(offsets * excess) / max(excess.sum(), 1)
    initial = [np.clip(background, 0, 255), 0, min(amplitude, 240),
               np.clip(center, -7.5, 7.5), 4.0, 1.5]
    solution = least_squares(
        lambda p: (_blurred_bar(p, offsets) - values) / noise, initial,
        bounds=([0, -2, 0, -9, 0.65, 0.45], [255, 2, 255, 9, 11, 5]),
        loss="soft_l1", f_scale=2.0, max_nfev=65,
    )
    parameters = solution.x
    rmse = np.sqrt(np.mean((_blurred_bar(parameters, offsets) - values) ** 2))
    return point + parameters[3] * normal, parameters, rmse


def _splines(distance, points, residual_budget):
    return [
        UnivariateSpline(distance, coordinate, s=len(distance) * residual_budget**2, k=3)
        for coordinate in points.T
    ]


def _at(splines, distance):
    return np.column_stack([spline(distance) for spline in splines])


def _tangents(splines, distance):
    directions = np.column_stack([spline.derivative()(distance) for spline in splines])
    norm = np.linalg.norm(directions, axis=1)
    if np.any(~np.isfinite(norm)) or np.any(norm <= 1e-12):
        raise ValueError("Undefined tangent in the fitted root centerline")
    return directions / norm[:, None]


def _fit_centers(gray, seed):
    distance = np.r_[0.0, np.cumsum(np.linalg.norm(np.diff(seed, axis=0), axis=1))]
    distance -= distance[-1]
    seed_splines = _splines(distance, seed, 0.45)
    samples = np.arange(-min(85.0, -distance[0] + 3), 12.001, 1.5)
    body = (samples >= -40) & (samples <= -15)
    if not body.any():
        raise ValueError("Seed path is too short for the transverse profile fit")
    inside = np.minimum(samples, -2.0)
    tangents = _tangents(seed_splines, inside)
    positions = _at(seed_splines, inside) + np.maximum(samples + 2.0, 0)[:, None] * tangents
    first_pass = [_fit_profile(gray, point, tangent) for point, tangent in zip(positions, tangents)]
    centers, parameters, errors = (np.asarray(items) for items in zip(*first_pass))
    peaks = parameters[:, 2] * erf(parameters[:, 4] / (np.sqrt(2) * parameters[:, 5]))
    widths = np.sqrt(parameters[:, 4]**2 / 3 + parameters[:, 5]**2)
    body_contrast = float(np.median(peaks[body]))
    body_width = float(np.median(widths[body]))
    valid = (peaks > max(5.0, body_contrast * 0.32)) & (widths > max(1.0, body_width * 0.35))
    valid &= (np.abs(parameters[:, 3]) < 8.8) & (errors < np.maximum(9.0, peaks * 0.45))
    # Two consecutive unsupported distal profiles end the candidate root body.
    breaks = np.flatnonzero((samples[:-1] > -10) & ~valid[:-1] & ~valid[1:])
    if breaks.size:
        valid[breaks[0]:] = False
    if valid.sum() < 20:
        raise ValueError("Insufficient original-intensity support for the root centerline")
    samples = samples[valid]
    centers = centers[valid]
    splines = _splines(samples, centers, 0.30)
    for index, (point, tangent) in enumerate(zip(_at(splines, samples), _tangents(splines, samples))):
        corrected, parameters, _ = _fit_profile(gray, point, tangent)
        centers[index] = corrected if abs(parameters[3]) < 4 and parameters[2] > 5 else point
    splines = _splines(samples, centers, 0.30)
    residual = np.linalg.norm(_at(splines, samples) - centers, axis=1)
    diagnostics = {
        "body_contrast": body_contrast, "body_profile_sigma_px": body_width,
        "fit_center_count": int(len(samples)), "profile_rmse_median": float(np.median(errors[valid])),
        "centerline_residual_rms_px": float(np.sqrt(np.mean(residual**2))),
    }
    return splines, samples, centers, diagnostics


def _optical_tip(gray, splines, samples):
    last = samples[-1]
    end = _at(splines, [last])[0]
    tangent = _tangents(splines, [last])[0]

    def curve(points):
        # Use a tangent continuation beyond the last illuminated transverse center.
        return _at(splines, np.minimum(points, last)) + np.maximum(points-last, 0)[:, None] * tangent

    grid = np.arange(-18, 20.001, 0.1)
    signal = _intensity(gray, curve(grid))
    derivative = np.gradient(gaussian_filter1d(signal, 10.0), 0.1)
    search = np.flatnonzero((grid >= -9) & (grid <= 14))
    index = search[np.argmin(derivative[search])]
    if not np.isfinite(derivative[index]) or derivative[index] >= 0:
        raise ValueError("No descending intensity edge at the root cap")
    left, middle, right = derivative[index-1:index+2]
    denominator = left - 2 * middle + right
    offset = 0.5 * (left - right) / denominator if abs(denominator) > 1e-12 else 0.0
    tip_s = grid[index] + np.clip(offset, -1, 1) * 0.1
    grid = np.r_[np.arange(max(-75, samples[0]), tip_s, 0.05), tip_s]
    return curve(grid), float(tip_s), float(derivative[index])


def _pixel_sizes(dpi):
    if dpi is None:
        raise ValueError("Source X/Y resolution in dots per inch is required")
    resolution = np.asarray(dpi, dtype=float)
    if resolution.ndim == 0:
        resolution = np.repeat(resolution, 2)
    if resolution.shape != (2,) or not np.isfinite(resolution).all() or np.any(resolution <= 0):
        raise ValueError("Source resolution must contain finite positive X/Y DPI values")
    return 25.4 / resolution


def fit_tip(gray, seed_xy, length_mm=DEFAULT_LENGTH_MM, dpi=None, mirror_x=False):
    """Return a source-patch centerline, tip, upstream point and signed chord angle.

    Supply a 2D 8-bit luminance patch from the original scan and a subpixel XY
    seed ordered from proximal body to distal cap in that same patch. The caller
    chooses mirror_x from the original plate orientation. Positive Y is down.
    Native-pixel profile widths and smoothing were tested on 1200-dpi scans.
    Fit success yields a candidate for review. ValueError records a failed fit.
    """
    pixel_mm = _pixel_sizes(dpi)
    length_mm = float(length_mm)
    if not np.isfinite(length_mm) or length_mm <= 0:
        raise ValueError("Root-tip segment length must be finite and positive")
    gray = np.asarray(gray)
    if gray.ndim != 2 or gray.dtype != np.uint8 or not gray.size:
        raise ValueError("Expected an 8-bit grayscale patch from the original scan")
    gray = gray.astype(float)
    seed = np.asarray(seed_xy, dtype=float)
    if seed.ndim != 2 or seed.shape[1] != 2 or len(seed) < 4 or not np.isfinite(seed).all():
        raise ValueError("Seed requires at least four finite original-image XY points")
    if np.any(np.linalg.norm(np.diff(seed, axis=0), axis=1) <= 0):
        raise ValueError("Seed path has repeated neighboring coordinates")
    bounds = np.array([gray.shape[1] - 1, gray.shape[0] - 1])
    if np.any(seed < 0) or np.any(seed > bounds):
        raise ValueError("Seed coordinates extend beyond the original image patch")

    splines, samples, centers, diagnostics = _fit_centers(gray, seed)
    curve, tip_s, gradient = _optical_tip(gray, splines, samples)
    steps = np.linalg.norm(np.diff(curve, axis=0) * pixel_mm, axis=1)
    if not len(steps) or not np.isfinite(steps).all() or np.any(steps <= 0):
        raise ValueError("Fitted curve has an undefined arc length")
    arc = np.r_[0.0, np.cumsum(steps)]
    back = arc[-1] - arc
    if length_mm > back[0]:
        raise ValueError(f"Requested {length_mm:g} mm exceeds fitted arc of {back[0]:.4g} mm")
    upstream = np.array([np.interp(length_mm, back[::-1], v[::-1]) for v in curve.T])
    tip = curve[-1]
    segment = np.vstack((upstream, curve[back <= length_mm]))
    if np.any(segment < 0) or np.any(segment > bounds):
        raise ValueError("Requested root-tip segment leaves the original image patch")
    dx, dy = (tip - upstream) * pixel_mm
    angle = float(np.degrees(np.arctan2(dy, -dx if mirror_x else dx)))
    diagnostics.update({
        "tip_s": tip_s, "last_fitted_center_s": float(samples[-1]),
        "tip_extrapolation_px": float(max(0.0, tip_s - samples[-1])),
        "tip_gradient": gradient, "tip_shift_from_mask_px": float(np.linalg.norm(tip - seed[-1])),
    })
    return {
        "angle_deg": angle, "tip_xy": tip, "upstream_xy": upstream,
        "curve_xy": curve, "back_arc_mm": back, "centers_xy": centers,
        "length_mm": length_mm, "diagnostics": diagnostics,
    }
