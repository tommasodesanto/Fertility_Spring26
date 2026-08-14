"use client";

import { useMemo, useState } from "react";
import liveStatus from "./generated-status.json";
import {
  AREAS,
  EDGES,
  NODES,
  OBJECTS,
  STATUS_COPY,
  TRACES,
  type MapEdge,
  type MapObject,
  type ObjectStatus,
} from "./project-data";

const VIEW_WIDTH = 1200;
const VIEW_HEIGHT = 760;
const STATUSES = Object.keys(STATUS_COPY) as ObjectStatus[];

type TargetRow = {
  name: string;
  target: number;
  standardError: number;
  weight: number;
  weightKind: string;
  status: string;
};

type ParameterRow = {
  name: string;
  lower?: number;
  upper?: number;
  value?: number;
  transform?: string;
  restriction: string;
};

const statusData = liveStatus as {
  generatedFrom: {
    repoRoot: string;
    calibrationStatusUpdated: string;
    targetSet: string;
    profileName: string;
    fingerprint: string;
  };
  current: {
    maintainedStrand: string;
    contractStatus: ObjectStatus;
    roomsCodeHold: boolean;
    closure: string;
    normalizedPrice: number | null;
    renewalPrice: number | null;
    stationaryScale: number | null;
    outsideFlow: number | null;
    reproductionRatio: number | null;
  };
  warnings: string[];
  targetRows: TargetRow[];
  parameters: ParameterRow[];
  sourceChecks: { path: string; exists: boolean }[];
};

const coords = Object.fromEntries(
  NODES.map((node) => [
    node.id,
    { x: ((node.x ?? 0) / 100) * VIEW_WIDTH, y: ((node.y ?? 0) / 100) * VIEW_HEIGHT },
  ]),
);

function edgePath(edge: MapEdge) {
  const start = coords[edge.from];
  const end = coords[edge.to];
  const midX = (start.x + end.x) / 2;
  const midY = (start.y + end.y) / 2;
  const dx = end.x - start.x;
  const dy = end.y - start.y;
  const length = Math.max(Math.hypot(dx, dy), 1);
  const bend = edge.bend ?? 0;
  const controlX = midX + (-dy / length) * bend;
  const controlY = midY + (dx / length) * bend;
  return `M ${start.x} ${start.y} Q ${controlX} ${controlY} ${end.x} ${end.y}`;
}

function statusSlug(status: ObjectStatus) {
  return status.replaceAll(" ", "-");
}

function humanize(value: string) {
  return value.replaceAll("_", " ").replace(/\b\w/g, (letter) => letter.toUpperCase());
}

function compactNumber(value: number | null | undefined, digits = 4) {
  if (value === null || value === undefined || !Number.isFinite(value)) return "not parsed";
  if (Math.abs(value) >= 1000) return value.toLocaleString("en-US", { maximumFractionDigits: 1 });
  return value.toLocaleString("en-US", { maximumFractionDigits: digits });
}

function sourceExists(path: string) {
  return statusData.sourceChecks.find((row) => row.path === path)?.exists;
}

export default function Home() {
  const [selectedId, setSelectedId] = useState("target-contract");
  const [traceId, setTraceId] = useState("calibration");
  const [paused, setPaused] = useState(false);
  const [statusFilter, setStatusFilter] = useState<ObjectStatus | null>(null);
  const [copied, setCopied] = useState<string | null>(null);

  const trace = TRACES.find((item) => item.id === traceId) ?? TRACES[0];
  const selected = OBJECTS[selectedId] ?? OBJECTS["target-contract"];
  const activeNodes = useMemo(() => new Set(trace.nodes), [trace.nodes]);
  const activeEdges = useMemo(() => new Set(trace.edges), [trace.edges]);

  async function copyText(key: string, text: string) {
    let succeeded = false;
    try {
      await navigator.clipboard.writeText(text);
      succeeded = true;
    } catch {
      const fallback = document.createElement("textarea");
      fallback.value = text;
      fallback.setAttribute("readonly", "");
      fallback.style.position = "fixed";
      fallback.style.opacity = "0";
      document.body.appendChild(fallback);
      fallback.select();
      succeeded = document.execCommand("copy");
      fallback.remove();
    }
    if (succeeded) {
      setCopied(key);
      window.setTimeout(() => setCopied((current) => (current === key ? null : current)), 1800);
    }
  }

  function selectTrace(nextTraceId: string) {
    const next = TRACES.find((item) => item.id === nextTraceId) ?? TRACES[0];
    setTraceId(next.id);
    setSelectedId(next.nodes[0]);
  }

  const prompt = [
    `Investigate the project object: ${selected.title}.`,
    `Economic meaning: ${selected.meaning}`,
    `Key equation: ${selected.equation}`,
    `Current trace: ${trace.label} — ${trace.question}`,
    `Start from these exact files:\n${selected.sources.map((path) => `- ${path}`).join("\n")}`,
    `Known caveats:\n${selected.caveats.map((item) => `- ${item}`).join("\n")}`,
    "Read AGENTS.md and the mandatory startup files first. Verify the implementation against the equation, trace inputs and outputs, and report any mismatch, weak identification, stale object, or missing reproducibility step. Do not edit source code unless I explicitly ask.",
  ].join("\n\n");

  return (
    <main className="shell">
      <header className="topbar">
        <div className="title-lockup">
          <p className="kicker">FERTILITY × HOUSING / RESEARCH SYSTEM</p>
          <h1>The economic argument, made inspectable.</h1>
          <p className="deck">
            Evidence, equilibrium, and demographic renewal shown as a small set of real contracts—not a graph of every file.
          </p>
        </div>
        <div className="live-card">
          <div className="live-line">
            <span className="pulse" aria-hidden="true" />
            <strong>Maintained strand</strong>
          </div>
          <span>{statusData.current.maintainedStrand}</span>
          <span className={`status-pill ${statusSlug(statusData.current.contractStatus)}`}>
            {statusData.current.contractStatus}
          </span>
          <small>Status refreshed from CALIBRATION_STATUS.md · {statusData.generatedFrom.calibrationStatusUpdated}</small>
        </div>
      </header>

      <section className="control-deck" aria-label="Map controls">
        <div className="trace-controls">
          <span className="control-label">TRACE MODE</span>
          {TRACES.map((item) => (
            <button
              key={item.id}
              className={trace.id === item.id ? "active" : ""}
              onClick={() => selectTrace(item.id)}
              aria-pressed={trace.id === item.id}
            >
              <span>{item.shortLabel}</span>
              {item.label}
            </button>
          ))}
        </div>
        <button className="motion-toggle" onClick={() => setPaused((value) => !value)} aria-pressed={paused}>
          <span>{paused ? "▶" : "Ⅱ"}</span> {paused ? "Resume packets" : "Pause packets"}
        </button>
      </section>

      <section className="trace-brief">
        <div>
          <span className="trace-code">{trace.shortLabel} / {String(TRACES.indexOf(trace) + 1).padStart(2, "0")}</span>
          <p>{trace.description}</p>
        </div>
        <blockquote>{trace.question}</blockquote>
      </section>

      <section className="workspace">
        <div className="map-wrap">
          <div className="map-meta">
            <div>
              <span>ACTIVE SYSTEM MAP</span>
              <small>{trace.nodes.length} blocks · {trace.packets.length} packet types in focus</small>
            </div>
            <div className="legend" aria-label="Status legend">
              {STATUSES.map((status) => (
                <button
                  key={status}
                  className={`${statusSlug(status)} ${statusFilter === status ? "selected" : ""}`}
                  title={STATUS_COPY[status]}
                  onClick={() => setStatusFilter((current) => (current === status ? null : status))}
                >
                  <span /> {status}
                </button>
              ))}
            </div>
          </div>

          <div className="mobile-flow" aria-label={`${trace.label} mobile trace`}>
            <div className={`mobile-packet-track ${paused ? "paused" : ""}`}><span /></div>
            <ol>
              {trace.nodes.map((nodeId, index) => {
                const node = OBJECTS[nodeId];
                return (
                  <li key={node.id}>
                    <button onClick={() => setSelectedId(node.id)} className={selectedId === node.id ? "selected" : ""}>
                      <span>{String(index + 1).padStart(2, "0")}</span>
                      <strong>{node.title}</strong>
                      <em className={statusSlug(node.status)}>{node.status}</em>
                    </button>
                  </li>
                );
              })}
            </ol>
          </div>

          <div className="map-scroll">
            <div className={`map-stage ${paused ? "paused" : ""}`}>
              {AREAS.map((area) => (
                <div
                  key={area.label}
                  className="area-band"
                  style={{ left: `${area.x}%`, width: `${area.width}%` }}
                >
                  <span>{area.label}</span>
                </div>
              ))}

              <svg className="edge-layer" viewBox={`0 0 ${VIEW_WIDTH} ${VIEW_HEIGHT}`} aria-hidden="true">
                <defs>
                  <marker id="arrow" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto">
                    <path d="M 0 0 L 7 3.5 L 0 7 z" />
                  </marker>
                  <marker id="arrow-active" markerWidth="7" markerHeight="7" refX="6" refY="3.5" orient="auto">
                    <path d="M 0 0 L 7 3.5 L 0 7 z" />
                  </marker>
                </defs>
                {EDGES.map((edge) => {
                  const isActive = activeEdges.has(edge.id);
                  return (
                    <path
                      key={edge.id}
                      id={edge.id}
                      d={edgePath(edge)}
                      className={`edge ${isActive ? "active" : ""}`}
                      markerEnd={`url(#${isActive ? "arrow-active" : "arrow"})`}
                    >
                      <title>{edge.label}</title>
                    </path>
                  );
                })}
                {trace.packets.map((packetId, index) => {
                  const packet = OBJECTS[packetId];
                  const edge = EDGES.find((item) => item.id === trace.edges[index % trace.edges.length]) ?? EDGES[0];
                  const start = coords[edge.from];
                  if (paused) {
                    return (
                      <circle
                        key={`${packetId}-${index}`}
                        cx={start.x}
                        cy={start.y}
                        r="7"
                        className={`moving-packet ${statusSlug(packet.status)}`}
                        onClick={() => setSelectedId(packetId)}
                      ><title>{packet.title}</title></circle>
                    );
                  }
                  return (
                    <circle
                      key={`${packetId}-${index}`}
                      r="7"
                      tabIndex={0}
                      role="button"
                      className={`moving-packet ${statusSlug(packet.status)}`}
                      onClick={() => setSelectedId(packetId)}
                      onKeyDown={(event) => {
                        if (event.key === "Enter" || event.key === " ") setSelectedId(packetId);
                      }}
                    >
                      <title>{packet.title}</title>
                      <animateMotion
                        path={edgePath(edge)}
                        dur={`${5.5 + (index % 4) * 1.15}s`}
                        begin={`${-(index * 0.85)}s`}
                        repeatCount="indefinite"
                      />
                    </circle>
                  );
                })}
              </svg>

              {NODES.map((node) => {
                const inTrace = activeNodes.has(node.id);
                const filtered = statusFilter !== null && node.status !== statusFilter;
                return (
                  <button
                    key={node.id}
                    className={`map-node ${selectedId === node.id ? "selected" : ""} ${inTrace ? "in-trace" : "off-trace"} ${filtered ? "filtered" : ""}`}
                    style={{ left: `${node.x}%`, top: `${node.y}%` }}
                    onClick={() => setSelectedId(node.id)}
                    aria-label={`${node.title}, ${node.status}`}
                  >
                    <span className="node-eyebrow">{node.eyebrow}</span>
                    <strong>{node.title}</strong>
                    <span className={`node-status ${statusSlug(node.status)}`}><i />{node.status}</span>
                  </button>
                );
              })}
            </div>
          </div>

          <div className="packet-key">
            <span>PACKETS IN THIS TRACE</span>
            <div>
              {trace.packets.map((packetId) => {
                const packet = OBJECTS[packetId];
                return (
                  <button key={packetId} onClick={() => setSelectedId(packetId)} className={selectedId === packetId ? "selected" : ""}>
                    <i className={statusSlug(packet.status)} /> {packet.title}
                  </button>
                );
              })}
            </div>
          </div>
        </div>

        <Inspector
          object={selected}
          traceLabel={trace.label}
          copied={copied}
          copyText={copyText}
          prompt={prompt}
        />
      </section>

      <footer className="footer-note">
        <span>MAP CONTRACT</span>
        <p>
          Status, target values, weights, bounds, fingerprints, and closure diagnostics are refreshed from canonical repository files. Economic descriptions and graph topology are intentionally curated.
        </p>
        <code>{statusData.generatedFrom.targetSet} · {statusData.generatedFrom.fingerprint.slice(0, 12)}…</code>
      </footer>
    </main>
  );
}

function Inspector({
  object,
  traceLabel,
  copied,
  copyText,
  prompt,
}: {
  object: MapObject;
  traceLabel: string;
  copied: string | null;
  copyText: (key: string, text: string) => Promise<void>;
  prompt: string;
}) {
  const [promptOpen, setPromptOpen] = useState(false);

  return (
    <aside className="inspector" aria-live="polite">
      <div className="inspector-topline">
        <span>{object.kind.toUpperCase()} / {object.area.toUpperCase()}</span>
        <span className={`status-pill ${statusSlug(object.status)}`}>{object.status}</span>
      </div>
      <h2>{object.title}</h2>
      <p className="meaning">{object.meaning}</p>

      <section className="detail-block equation-block">
        <h3>KEY EQUATION</h3>
        <code>{object.equation}</code>
      </section>

      <section className="io-grid">
        <div>
          <h3>INPUTS</h3>
          {object.inputs.map((item) => <span key={item}>{item}</span>)}
        </div>
        <div>
          <h3>OUTPUTS</h3>
          {object.outputs.map((item) => <span key={item}>{item}</span>)}
        </div>
      </section>

      {object.id === "target-contract" && <TargetContract />}
      {object.id === "parameters" && <ParameterContract />}
      {object.id === "population-closure" && <ClosureSnapshot />}

      <section className="detail-block sources-block">
        <h3>EXACT SOURCE FILES</h3>
        <ul>
          {object.sources.map((path, index) => {
            const exists = sourceExists(path);
            return (
              <li key={path}>
                <span className={exists === false ? "missing" : ""}>
                  <i>{exists === false ? "!" : String(index + 1).padStart(2, "0")}</i>
                  <code>{path}</code>
                </span>
                <button onClick={() => copyText(`source-${index}`, path)} aria-label={`Copy ${path}`}>
                  {copied === `source-${index}` ? "copied" : "copy"}
                </button>
              </li>
            );
          })}
        </ul>
      </section>

      <section className="detail-block caveat-block">
        <h3>CURRENT STATUS + CAVEATS</h3>
        <p className={`status-explainer ${statusSlug(object.status)}`}>{STATUS_COPY[object.status]}</p>
        {object.caveats.map((item) => <p key={item}>{item}</p>)}
      </section>

      <section className="detail-block command-block">
        <h3>REPRODUCE / INSPECT</h3>
        <div>
          <code>{object.command}</code>
          <button onClick={() => copyText("command", object.command)}>
            {copied === "command" ? "Copied" : "Copy command"}
          </button>
        </div>
      </section>

      <button className="prompt-button" onClick={() => setPromptOpen((value) => !value)} aria-expanded={promptOpen}>
        <span>
          <small>AI HANDOFF / {traceLabel.toUpperCase()}</small>
          {promptOpen ? "Hide investigation prompt" : "Generate investigation prompt"}
        </span>
        <b>{promptOpen ? "−" : "↗"}</b>
      </button>
      {promptOpen && (
        <section className="prompt-output" aria-label="Generated AI investigation prompt">
          <code>{prompt}</code>
          <button onClick={() => copyText("prompt", prompt)}>
            {copied === "prompt" ? "Prompt copied" : "Copy prompt"}
          </button>
        </section>
      )}
    </aside>
  );
}

function TargetContract() {
  return (
    <details className="data-contract" open>
      <summary>
        <span>LIVE TARGET CONTRACT</span>
        <em>{statusData.targetRows.length} rows</em>
      </summary>
      <div className="table-scroll">
        <table>
          <thead><tr><th>Moment</th><th>Target</th><th>SE</th><th>Weight</th></tr></thead>
          <tbody>
            {statusData.targetRows.map((row) => (
              <tr key={row.name} className={row.status === "empirically held" ? "held-row" : ""}>
                <td title={row.weightKind}>{humanize(row.name)}</td>
                <td>{compactNumber(row.target, 5)}</td>
                <td>{compactNumber(row.standardError, 5)}</td>
                <td>{compactNumber(row.weight, 1)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <p>Target values and weights are parsed from <code>e5_profile.py</code>. This table is the contract, not a model-fit readout.</p>
    </details>
  );
}

function ParameterContract() {
  return (
    <details className="data-contract" open>
      <summary>
        <span>PARAMETER CONTRACT</span>
        <em>{statusData.parameters.length} rows</em>
      </summary>
      <div className="table-scroll">
        <table>
          <thead><tr><th>Parameter</th><th>Restriction</th><th>Range / value</th></tr></thead>
          <tbody>
            {statusData.parameters.map((row) => (
              <tr key={row.name}>
                <td>{row.name}</td>
                <td>{row.restriction}</td>
                <td>{row.restriction === "estimated" ? `[${compactNumber(row.lower)}, ${compactNumber(row.upper)}]` : compactNumber(row.value)}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <p>Bounds and external restrictions are parsed from <code>e5_profile.py</code>; current estimates are intentionally not shown without a complete fit packet.</p>
    </details>
  );
}

function ClosureSnapshot() {
  const rows = [
    ["Normalized price", statusData.current.normalizedPrice],
    ["Renewal price", statusData.current.renewalPrice],
    ["Stationary scale", statusData.current.stationaryScale],
    ["Outside flow M", statusData.current.outsideFlow],
    ["Effective B₀/E₀", statusData.current.reproductionRatio],
  ] as const;
  return (
    <section className="closure-snapshot">
      <div><span>PARSED STATIONARY DIAGNOSTIC</span><em>{statusData.generatedFrom.calibrationStatusUpdated}</em></div>
      <dl>
        {rows.map(([label, value]) => <div key={label}><dt>{label}</dt><dd>{compactNumber(value, 8)}</dd></div>)}
      </dl>
      <p>Read as an open stationary accounting result, not a transition forecast.</p>
    </section>
  );
}
