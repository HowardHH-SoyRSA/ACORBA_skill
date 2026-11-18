// Use the spreadsheet runtime bundled with Codex Desktop.
import fs from "node:fs/promises";
import path from "node:path";
import { pathToFileURL } from "node:url";

const [, , packageRoot, output] = process.argv;
const pkg = JSON.parse(await fs.readFile(path.join(packageRoot, "package.json"), "utf8"));
const entry = path.join(packageRoot, pkg.exports["."]);
const { Workbook, SpreadsheetFile } = await import(pathToFileURL(entry).href);
const chunks = [];
for await (const chunk of process.stdin) chunks.push(chunk);
const sheets = JSON.parse(Buffer.concat(chunks).toString("utf8"));
const workbook = Workbook.create();

for (const { name, data, header, widths, decimal_cols } of sheets) {
  const sheet = workbook.worksheets.add(name);
  const end = String.fromCharCode(64 + header.length);
  const last = data.length + 1;
  sheet.getRange(`A1:${end}${last}`).values = [header, ...data];
  const heading = sheet.getRange(`A1:${end}1`);
  heading.format.fill = "#254C58";
  heading.format.font.color = "#FFFFFF";
  heading.format.font.bold = true;
  heading.format.rowHeight = 24;
  for (const [index, width] of widths.entries()) {
    const col = String.fromCharCode(65 + index);
    sheet.getRange(`${col}1:${col}${last}`).format.columnWidth = width;
    if (decimal_cols.includes(index)) {
      sheet.getRange(`${col}2:${col}${last}`).setNumberFormat("0.00");
    }
  }
  sheet.freezePanes.freezeRows(1);
}
await (await SpreadsheetFile.exportXlsx(workbook)).save(output);
