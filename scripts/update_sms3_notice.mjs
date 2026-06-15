#!/usr/bin/env node
import { readdirSync, readFileSync, writeFileSync } from "node:fs";
import path from "node:path";
import { fileURLToPath } from "node:url";

const scriptDir = path.dirname(fileURLToPath(import.meta.url));
const repoRoot = path.resolve(scriptDir, "..");
const docsDir = path.join(repoRoot, "docs");
const stylesheetPath = path.join(docsDir, "styles", "stylesheet.css");
const sms3Url = "https://paulstothard.github.io/sequence-manipulation-suite/";

const noticeStart = "                <!-- SMS3 NOTICE START -->";
const noticeEnd = "                <!-- SMS3 NOTICE END -->";
const noticeBlock = [
  noticeStart,
  "                <tr>",
  '                  <td class="description">',
  '                    <div class="sms3_notice">',
  "                      <b>New version available:</b>",
  `                      <a href="${sms3Url}">Sequence Manipulation Suite 3 (SMS3)</a>`,
  "                    </div>",
  "                  </td>",
  "                </tr>",
  noticeEnd
].join("\n");

const noticePattern = new RegExp(
  `\\n${escapeRegExp(noticeStart)}[\\s\\S]*?${escapeRegExp(noticeEnd)}`,
  "g"
);
const titleRowPattern =
  /(\n                <tr>\n                  <td class="title">[\s\S]*?<\/td>\n                <\/tr>)/;

const cssStart = "/* SMS3 NOTICE START */";
const cssEnd = "/* SMS3 NOTICE END */";
const cssBlock = [
  cssStart,
  "div.sms3_notice {",
  "  margin: 0.75em 0 1em 0;",
  "  padding: 0.6em 0.75em;",
  "  border: 1px solid #cccccc;",
  "  color: #000000;",
  "  background-color: #f6f6f6;",
  "}",
  "div.sms3_notice a {",
  "  color: #000099;",
  "  text-decoration: none;",
  "}",
  "div.sms3_notice a:visited {",
  "  color: #000099;",
  "  text-decoration: none;",
  "}",
  "div.sms3_notice a:hover {",
  "  color: #ff0000;",
  "  text-decoration: underline;",
  "}",
  "div.sms3_notice a:active {",
  "  color: #000099;",
  "  text-decoration: none;",
  "}",
  cssEnd
].join("\n");
const cssPattern = new RegExp(
  `\\n?${escapeRegExp(cssStart)}[\\s\\S]*?${escapeRegExp(cssEnd)}\\n?`,
  "g"
);

function escapeRegExp(value) {
  return value.replace(/[.*+?^${}()|[\]\\]/g, "\\$&");
}

function updateHtmlFile(filePath) {
  const original = readFileSync(filePath, "utf8");
  const withoutExistingNotice = original.replace(noticePattern, "");
  if (!titleRowPattern.test(withoutExistingNotice)) {
    throw new Error(`Could not find title row in ${path.relative(repoRoot, filePath)}`);
  }
  const updated = withoutExistingNotice.replace(titleRowPattern, `$1\n${noticeBlock}`);
  if (updated !== original) {
    writeFileSync(filePath, updated, "utf8");
    return true;
  }
  return false;
}

function updateStylesheet() {
  const original = readFileSync(stylesheetPath, "utf8");
  const withoutExistingNotice = original.replace(cssPattern, "").replace(/\s+$/u, "");
  const updated = `${withoutExistingNotice}\n${cssBlock}\n`;
  if (updated !== original) {
    writeFileSync(stylesheetPath, updated, "utf8");
    return true;
  }
  return false;
}

const htmlFiles = readdirSync(docsDir)
  .filter((fileName) => fileName.endsWith(".html"))
  .sort()
  .map((fileName) => path.join(docsDir, fileName));

const changedHtmlCount = htmlFiles.filter(updateHtmlFile).length;
const changedCss = updateStylesheet();

console.log(
  `Updated SMS3 notice in ${changedHtmlCount} HTML file(s)` +
    (changedCss ? " and stylesheet." : ".")
);
