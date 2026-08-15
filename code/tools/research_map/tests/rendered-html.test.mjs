import assert from "node:assert/strict";
import { readdir, readFile } from "node:fs/promises";
import test from "node:test";

test("production bundle contains the research map", async () => {
  const dist = new URL("../dist/", import.meta.url);
  const files = (await readdir(dist, { recursive: true })).filter((path) => path.endsWith(".js"));
  const bundle = (await Promise.all(files.map((path) => readFile(new URL(path, dist), "utf8")))).join("\n");
  assert.ok(bundle.includes("Fertility × Housing Research Map"));
  assert.ok(bundle.includes("The economic argument, made inspectable."));
  assert.ok(bundle.includes("Five economic operations, one fixed point."));
  assert.ok(bundle.includes("IN PLAIN ENGLISH"));
  assert.ok(bundle.includes("LIVE RESEARCH CONTRACT"));
  assert.ok(bundle.includes("HOW THIS BLOCK WORKS"));
  assert.ok(bundle.includes("WORKED EXAMPLE"));
  assert.ok(bundle.includes("JUMP TO AN OBJECT"));
  assert.ok(bundle.includes("Housing-market equilibrium"));
  assert.doesNotMatch(bundle, /codex-preview|SkeletonPreview|Your site is taking shape/i);
});

test("generated status carries the live contract", async () => {
  const payload = JSON.parse(await readFile(new URL("../app/generated-status.json", import.meta.url), "utf8"));
  assert.notEqual(payload.generatedFrom.targetSet, "unknown");
  assert.notEqual(payload.generatedFrom.profileName, "unknown");
  assert.notEqual(payload.generatedFrom.fingerprint, "unknown");
  assert.match(payload.generatedFrom.repositoryUrl, /^https:\/\/github\.com\//);
  assert.notEqual(payload.generatedFrom.branch, "unknown");
  assert.notEqual(payload.generatedFrom.commit, "unknown");
  assert.equal(payload.targetRows.length, 12);
  assert.equal(payload.parameters.filter((row) => row.restriction === "estimated").length, 10);
  assert.equal(payload.current.roomsCodeHold, true);
  assert.equal(payload.current.closure, "minimal fixed-inflow renewal closure");
  assert.equal(payload.current.calendarTransitionSolved, false);
  assert.equal(payload.sourceChecks.find((row) => row.path === "CALIBRATION_STATUS.md")?.exists, true);
  assert.equal(payload.sourceChecks.find((row) => row.path.endsWith("/e5_profile.py"))?.exists, true);
});
