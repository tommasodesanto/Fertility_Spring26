#!/usr/bin/env node

import { existsSync, readFileSync, writeFileSync } from "node:fs";
import { execFileSync } from "node:child_process";
import { dirname, resolve } from "node:path";
import { fileURLToPath } from "node:url";

const toolDir = resolve(dirname(fileURLToPath(import.meta.url)), "..");
const repoFlag = process.argv.indexOf("--repo");
const repoRoot = repoFlag >= 0 ? resolve(process.argv[repoFlag + 1]) : resolve(toolDir, "../../..");

const read = (relative) => readFileSync(resolve(repoRoot, relative), "utf8");
const statusText = read("CALIBRATION_STATUS.md");
const e5Text = read("code/model/intergen_eqscale_seq_optimized/e5_profile.py");
const oldCalibrationText = read("code/model/intergen_eqscale_seq_optimized/calibration.py");
const launcherText = read("code/cluster/submit_intergen_e5b_idfe.sh");

function git(...args) {
  try {
    return execFileSync("git", args, { cwd: repoRoot, encoding: "utf8" }).trim();
  } catch {
    return "unknown";
  }
}

function normalizeRemote(remote) {
  if (remote === "unknown") return "unknown";
  if (remote.startsWith("git@github.com:")) {
    return `https://github.com/${remote.slice("git@github.com:".length).replace(/\.git$/, "")}`;
  }
  return remote.replace(/\.git$/, "");
}

function match(text, expression, fallback = null) {
  const result = text.match(expression);
  return result ? result[1] : fallback;
}

function block(text, name) {
  const startExpression = new RegExp(`^${name}[^=]*=\\s*\\{`, "m");
  const startMatch = startExpression.exec(text);
  if (!startMatch) throw new Error(`Could not locate ${name}`);
  const start = startMatch.index + startMatch[0].length;
  let depth = 1;
  let quote = null;
  for (let index = start; index < text.length; index += 1) {
    const character = text[index];
    const previous = text[index - 1];
    if ((character === '"' || character === "'") && previous !== "\\") {
      quote = quote === character ? null : quote === null ? character : quote;
    }
    if (quote) continue;
    if (character === "{") depth += 1;
    if (character === "}") depth -= 1;
    if (depth === 0) return text.slice(start, index);
  }
  throw new Error(`Unclosed ${name}`);
}

function numericDictionary(text, name, references = {}) {
  const result = {};
  for (const line of block(text, name).split("\n")) {
    const row = line.match(/^\s*["']([^"']+)["']\s*:\s*([^,#]+)(?:,|$)/);
    if (!row) continue;
    const key = row[1];
    const raw = row[2].trim();
    if (/^-?\d+(?:\.\d+)?(?:e[+-]?\d+)?$/i.test(raw)) {
      result[key] = Number(raw);
      continue;
    }
    const reference = raw.match(/^([A-Z0-9_]+)\[["']([^"']+)["']\]$/);
    if (reference && references[reference[1]]?.[reference[2]] !== undefined) {
      result[key] = references[reference[1]][reference[2]];
    }
  }
  return result;
}

function stringConstant(text, name, fallback) {
  return match(text, new RegExp(`^${name}\\s*=\\s*["']([^"']+)["']`, "m"), fallback);
}

const oldTargets = numericDictionary(oldCalibrationText, "OLD_NONLOCATION_TARGETS");
const targets = numericDictionary(e5Text, "E5_TARGETS", { OLD_NONLOCATION_TARGETS: oldTargets });
const measuredSes = numericDictionary(e5Text, "E5_MEASURED_SES");

const targetRows = Object.entries(targets).map(([name, target]) => {
  const measured = measuredSes[name] !== undefined;
  const standardError = measured ? measuredSes[name] : 0.05 * Math.abs(target);
  return {
    name,
    target,
    standardError,
    weight: 1 / standardError ** 2,
    weightKind: measured ? "measured / declared" : "synthetic 5%",
    status: name === "housing_increment_0to1" ? "empirically held" : "provisional",
  };
});

const domainText = match(e5Text, /E5_DOMAIN[^=]*=\s*\(([^]*?)\n\)/m, "");
const parameters = [];
for (const row of domainText.matchAll(/\(["']([^"']+)["'],\s*([-\d.]+),\s*([-\d.]+),\s*["']([^"']+)["']\)/g)) {
  parameters.push({ name: row[1], lower: Number(row[2]), upper: Number(row[3]), transform: row[4], restriction: "estimated" });
}

const fixedText = block(e5Text, "E5_FIXED");
for (const row of fixedText.matchAll(/["']([^"']+)["']\s*:\s*([-\d.]+)/g)) {
  parameters.push({ name: row[1], value: Number(row[2]), restriction: "externally fixed" });
}

const exactSourcePaths = [
  "CALIBRATION_STATUS.md",
  "code/cluster/submit_intergen_e5b_idfe.sh",
  "code/cluster/submit_intergen_e5b_idfe_collector.sh",
  "code/data/cps_fertility/build_cps_fertility_targets.py",
  "code/data/mms_center_periphery/build_intergen_one_market_housing_targets.R",
  "code/data/moment_standard_errors/SPEC_moment_bootstrap_se.md",
  "code/data/moment_standard_errors/build_moment_bootstrap_se.R",
  "code/data/nchs_natality_timing/build_first_birth_timing.R",
  "code/data/psid_followup_mar2026/sa_rooms_first_birth_one_variant_v1.do",
  "code/model/dt_cp_model/objective.py",
  "code/model/dt_cp_model/solver.py",
  "code/model/intergen_eqscale_seq_optimized/e5_profile.py",
  "code/model/intergen_eqscale_seq_optimized/solver.py",
  "code/model/intergen_eqscale_seq_optimized/calibration.py",
  "code/model/intergen_eqscale_seq_optimized/collect_e1.py",
  "code/model/intergen_eqscale_seq_optimized/externals.py",
  "code/model/intergen_eqscale_seq_optimized/parameters.py",
  "code/model/intergen_eqscale_seq_optimized/run_e1_chain.py",
  "code/model/intergen_eqscale_seq_optimized/target_system.py",
  "code/model/intergen_eqscale_seq_optimized/tests/test_bequest_target_moments.py",
  "code/model/intergen_eqscale_seq_optimized/tests/test_eqscale_seq.py",
  "code/model/tools/build_population_closure_update.py",
  "code/model/tools/audit_closed_reproductive_closure.py",
  "code/model/tools/build_intergen_mechanics_packet.py",
  "code/model/tools/check_population_closure.py",
  "code/model/tools/run_e5_repaired_policy_with_entry.py",
  "code/model/tools/run_intergen_funded_policy_with_entry.py",
  "latex/dynamic_intergenerational_housing_fertility_model.tex",
  "latex/model_writeup.tex",
  "output/model/population_closure_update/",
  "output/pdf/dynamic_intergenerational_housing_fertility_model.pdf",
];

const data = {
  generatedFrom: {
    repoRoot,
    calibrationStatusUpdated: match(statusText, /^Updated:\s*`([^`]+)`/m, "unknown"),
    targetSet: stringConstant(e5Text, "E5_TARGET_SET", "unknown"),
    profileName: stringConstant(e5Text, "E5_PROFILE_NAME", "unknown"),
    fingerprint: match(launcherText, /^export E5_EXPECTED_TARGET_FINGERPRINT=([^\s]+)/m, "unknown"),
    repositoryUrl: normalizeRemote(git("config", "--get", "remote.origin.url")),
    branch: git("branch", "--show-current"),
    commit: git("rev-parse", "HEAD"),
  },
  current: {
    maintainedStrand: match(
      statusText,
      /^The ([^.\n]+) has been rerun under/m,
      "maintained model strand not parsed",
    ),
    contractStatus: /remain[^\n]*under empirical hold|under empirical hold/i.test(statusText)
      ? "empirically held"
      : "provisional",
    roomsCodeHold: /ACTUALROOMS_[^]*?under empirical hold/i.test(statusText),
    closure: "open stationary renewal wrapper",
    normalizedPrice: Number(match(statusText, /normalized and\s+renewal asset prices are respectively `([\d.]+)`/m, "NaN")),
    renewalPrice: Number(match(statusText, /renewal asset prices are respectively `[\d.]+` and\s+`([\d.]+)`/m, "NaN")),
    stationaryScale: Number(match(statusText, /renewal scale is `([\d.]+)`/m, "NaN")),
    outsideFlow: Number(match(statusText, /fixed outside\s+flow `M=E0-B0=([\d.]+)`/m, "NaN")),
    reproductionRatio: Number(match(statusText, /effective reproduction is `B0\/E0 = ([\d.]+)`/m, "NaN")),
  },
  warnings: [
    "The first-birth rooms target is under an empirical valid-code hold; no refit or policy result using it is promoted.",
    "The maintained shared child clock differs from the proposed independent-child maturation architecture; changing it requires re-estimation.",
    "The renewal result is an open stationary comparison, not a calendar-time U.S. transition.",
    "Historical spatial and pre-correction losses are not comparable to the current one-market target system.",
  ],
  targetRows,
  parameters,
  sourceChecks: exactSourcePaths.map((path) => ({ path, exists: existsSync(resolve(repoRoot, path)) })),
};

const output = resolve(toolDir, "app/generated-status.json");
writeFileSync(output, `${JSON.stringify(data, null, 2)}\n`);
console.log(`Refreshed ${output}`);
console.log(`Target set: ${data.generatedFrom.targetSet}; ${targetRows.length} targets; ${parameters.length} parameter rows.`);
