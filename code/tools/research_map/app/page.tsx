"use client";

import katex from "katex";
import { useEffect, useMemo, useState } from "react";
import { EXPLANATIONS } from "./explanations";
import liveStatus from "./generated-status.json";
import { RICH_CONTENT } from "./rich-content";
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
  type Trace,
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
    repositoryUrl: string;
    branch: string;
    commit: string;
  };
  current: {
    maintainedStrand: string;
    contractStatus: ObjectStatus;
    roomsCodeHold: boolean;
    closure: string;
    calendarTransitionSolved: boolean;
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

function repositoryLink(path: string) {
  if (statusData.generatedFrom.repositoryUrl === "unknown") return null;
  const revision = statusData.generatedFrom.commit === "unknown" ? statusData.generatedFrom.branch : statusData.generatedFrom.commit;
  const mode = path.endsWith("/") ? "tree" : "blob";
  return `${statusData.generatedFrom.repositoryUrl}/${mode}/${revision}/${path.replace(/^\//, "")}`;
}

function sourceKind(path: string) {
  if (path === "CALIBRATION_STATUS.md") return "canonical status";
  if (path.startsWith("code/data/")) return "empirical builder";
  if (path.includes("/tests/")) return "verification test";
  if (path.startsWith("code/cluster/")) return "cluster contract";
  if (path.startsWith("code/model/tools/")) return "diagnostic driver";
  if (path.startsWith("code/model/")) return "model source";
  if (path.startsWith("latex/")) return path.endsWith(".pdf") ? "rendered paper" : "paper source";
  if (path.startsWith("output/")) return "generated artifact";
  return "project source";
}

function MathExpression({ tex, display = false }: { tex: string; display?: boolean }) {
  const html = katex.renderToString(tex, {
    displayMode: display,
    throwOnError: false,
    strict: false,
    output: "htmlAndMathml",
  });
  return <span className={display ? "math-display" : "math-inline"} dangerouslySetInnerHTML={{ __html: html }} />;
}

const SYSTEM_STAGES = [
  {
    number: "01",
    title: "Measure",
    objectId: "empirical-builders",
    packet: "moments + uncertainty",
    description: "Survey estimators define the empirical objects the model must explain.",
  },
  {
    number: "02",
    title: "Choose",
    objectId: "household-problem",
    packet: "policy rules",
    description: "Households choose fertility, tenure, rooms, consumption, and saving over the lifecycle.",
  },
  {
    number: "03",
    title: "Aggregate",
    objectId: "kfe",
    packet: "distribution g(x)",
    description: "The forward equation turns individual rules into a stationary population distribution.",
  },
  {
    number: "04",
    title: "Clear",
    objectId: "price-clearing",
    packet: "demand ↔ price",
    description: "Physical housing demand meets supply; prices return to household budgets.",
  },
  {
    number: "05",
    title: "Renew",
    objectId: "population-closure",
    packet: "births + entrants",
    description: "Births mature into potential entrants and population accounting determines stationary scale.",
  },
];

function SystemGuide({ selectObject }: { selectObject: (id: string) => void }) {
  return (
    <section className="system-guide" aria-labelledby="system-guide-title">
      <div className="guide-heading">
        <span>START HERE / THE MODEL IN ONE LOOP</span>
        <h2 id="system-guide-title">Five economic operations, one fixed point.</h2>
        <p>
          Read left to right once to understand the mechanism. Then notice the return arrows: prices change choices, choices change the distribution, and births change future entry.
        </p>
      </div>
      <div className="guide-stages">
        {SYSTEM_STAGES.map((stage, index) => (
          <button key={stage.objectId} onClick={() => selectObject(stage.objectId)}>
            <span>{stage.number}</span>
            <strong>{stage.title}</strong>
            <p>{stage.description}</p>
            <em>{stage.packet}</em>
            {index < SYSTEM_STAGES.length - 1 && <i aria-hidden="true">→</i>}
          </button>
        ))}
      </div>
      <div className="fixed-point-note">
        <span>THE CLOSURE</span>
        <p><strong>This is not a one-way pipeline.</strong> A solution is reached only when household rules, the stationary distribution, housing prices, and renewal accounting are mutually consistent.</p>
      </div>
    </section>
  );
}

function TraceRoute({
  trace,
  selectedId,
  selectObject,
}: {
  trace: Trace;
  selectedId: string;
  selectObject: (id: string) => void;
}) {
  return (
    <section className="trace-route" aria-label={`${trace.label} causal sequence`}>
      <span>FOLLOW THIS TRACE</span>
      <div>
        {trace.nodes.map((nodeId, index) => {
          const node = OBJECTS[nodeId];
          return (
            <button key={node.id} className={selectedId === node.id ? "selected" : ""} onClick={() => selectObject(node.id)}>
              <i>{String(index + 1).padStart(2, "0")}</i>
              <strong>{node.title}</strong>
            </button>
          );
        })}
      </div>
    </section>
  );
}

function MechanismGraphic({ object }: { object: MapObject }) {
  return (
    <section className="mechanism-card" aria-label={`${object.title} input and output diagram`}>
      <div className="mechanism-side">
        <span>RECEIVES</span>
        {object.inputs.slice(0, 4).map((item) => <small key={item}>{item}</small>)}
      </div>
      <i aria-hidden="true">→</i>
      <div className="mechanism-core">
        <span>{object.kind}</span>
        <strong>{object.title}</strong>
        <small>{object.area}</small>
      </div>
      <i aria-hidden="true">→</i>
      <div className="mechanism-side output">
        <span>PRODUCES</span>
        {object.outputs.slice(0, 4).map((item) => <small key={item}>{item}</small>)}
      </div>
    </section>
  );
}

function ResearchSnapshot({ selectObject }: { selectObject: (id: string) => void }) {
  const measured = statusData.targetRows.filter((row) => row.weightKind === "measured / declared").length;
  const synthetic = statusData.targetRows.length - measured;
  const estimated = statusData.parameters.filter((row) => row.restriction === "estimated").length;
  const externallyFixed = statusData.parameters.length - estimated;
  const held = statusData.targetRows.filter((row) => row.status === "empirically held").length;

  const cards = [
    {
      id: "target-contract",
      number: String(statusData.targetRows.length).padStart(2, "0"),
      label: "Target rows",
      detail: `${measured} measured/declared SE · ${synthetic} synthetic · ${held} held`,
      kind: "observed",
    },
    {
      id: "parameters",
      number: String(estimated).padStart(2, "0"),
      label: "Free parameters",
      detail: `${externallyFixed} externally fixed · identification must be explicit`,
      kind: "estimated",
    },
    {
      id: "price-clearing",
      number: "p,g,S",
      label: "Joint equilibrium",
      detail: "household rules · distribution · price · stationary scale",
      kind: "endogenous",
    },
    {
      id: "population-closure",
      number: "OPEN",
      label: "Stationary boundary",
      detail: statusData.current.calendarTransitionSolved
        ? `${statusData.current.closure} · transition path recorded`
        : `${statusData.current.closure} · no calendar-time transition solved`,
      kind: "diagnostic",
    },
  ];

  return (
    <section className="research-snapshot" aria-labelledby="snapshot-title">
      <div className="snapshot-heading">
        <span>LIVE RESEARCH CONTRACT</span>
        <h2 id="snapshot-title">What is disciplined, chosen, cleared, and still open.</h2>
        <p>Counts come from the active target and parameter profile. Scope labels come from the canonical status contract.</p>
      </div>
      <div className="snapshot-grid">
        {cards.map((card) => (
          <button key={card.id} onClick={() => selectObject(card.id)}>
            <span className={`snapshot-number ${card.kind}`}>{card.number}</span>
            <strong>{card.label}</strong>
            <small>{card.detail}</small>
            <i aria-hidden="true">↗</i>
          </button>
        ))}
      </div>
    </section>
  );
}

function ObjectFinder({
  query,
  setQuery,
  selectObject,
}: {
  query: string;
  setQuery: (value: string) => void;
  selectObject: (id: string) => void;
}) {
  const results = useMemo(() => {
    const normalized = query.trim().toLowerCase();
    if (!normalized) return [];
    return [...NODES, ...Object.values(OBJECTS).filter((object) => object.kind === "packet")]
      .filter((object, index, array) => array.findIndex((candidate) => candidate.id === object.id) === index)
      .filter((object) => `${object.title} ${object.meaning} ${object.area} ${object.inputs.join(" ")} ${object.outputs.join(" ")}`.toLowerCase().includes(normalized))
      .slice(0, 7);
  }, [query]);

  function choose(id: string) {
    selectObject(id);
    setQuery("");
  }

  return (
    <section className="object-finder" aria-label="Find a research object">
      <label htmlFor="object-search">JUMP TO AN OBJECT</label>
      <div className="finder-input">
        <span aria-hidden="true">⌕</span>
        <input
          id="object-search"
          value={query}
          onChange={(event) => setQuery(event.target.value)}
          placeholder="Try ‘bequest’, ‘rooms’, ‘price’, or ‘entrant’"
          autoComplete="off"
        />
        {query && <button onClick={() => setQuery("")} aria-label="Clear object search">×</button>}
      </div>
      {query && (
        <div className="finder-results">
          {results.length > 0 ? results.map((object) => (
            <button key={object.id} onClick={() => choose(object.id)}>
              <span>{object.area}</span>
              <strong>{object.title}</strong>
              <small className={statusSlug(object.status)}>{object.status}</small>
            </button>
          )) : <p>No matching research object.</p>}
        </div>
      )}
    </section>
  );
}

function ContextGraphic({ object }: { object: MapObject }) {
  const caption = <div className="graphic-caption"><span>ECONOMIC VIEW</span><small>conceptual · not to scale</small></div>;

  if (object.area === "Evidence") {
    return (
      <section className="context-graphic evidence-graphic">
        {caption}
        <div className="pipeline-graphic">
          <div><span>01</span><strong>Records</strong><small>people · households · births</small></div><i>→</i>
          <div><span>02</span><strong>Estimator</strong><small>sample · FE · weights · cluster</small></div><i>→</i>
          <div><span>03</span><strong>Target object</strong><small>estimate · SE · provenance</small></div>
        </div>
      </section>
    );
  }

  if (object.area === "Calibration") {
    return (
      <section className="context-graphic calibration-graphic">
        {caption}
        <div className="loop-graphic">
          <div><b>θ</b><small>candidate</small></div><i>→</i>
          <div><b>Solve</b><small>p, policies, g</small></div><i>→</i>
          <div><b>m(θ)</b><small>model moments</small></div><i>→</i>
          <div><b>J(θ)</b><small>row gaps</small></div>
          <span>search updates θ only after a complete equilibrium</span>
        </div>
      </section>
    );
  }

  if (object.area === "Households") {
    return (
      <section className="context-graphic household-graphic">
        {caption}
        <div className="state-vector"><span>age</span><span>assets</span><span>income</span><span>tenure</span><span>parity</span><span>children home</span></div>
        <div className="lifecycle-line">
          <div><b>Enter</b><small>inherit · work</small></div>
          <div><b>Family years</b><small>birth · space · tenure</small></div>
          <div><b>Later life</b><small>save · own · dissave</small></div>
          <div><b>Exit</b><small>estate · bequest</small></div>
        </div>
      </section>
    );
  }

  if (object.area === "Aggregation") {
    return (
      <section className="context-graphic distribution-graphic">
        {caption}
        <div className="mass-flow">
          <div className="mass-cloud"><i /><i /><i /><i /><i /><i /><i /><i /></div>
          <span>policies + transitions</span>
          <div className="mass-bars"><i /><i /><i /><i /><i /></div>
        </div>
        <p>Individual decisions become a distribution only after every state transition is applied forward.</p>
      </section>
    );
  }

  if (object.area === "Markets") {
    return (
      <section className="context-graphic market-graphic">
        {caption}
        <div className="market-chart" aria-label="Conceptual housing supply and demand crossing">
          <span className="axis-y">rent / user cost</span><span className="axis-x">physical housing</span>
          <i className="demand-line" /><i className="supply-line" /><b className="market-dot" />
          <small className="demand-label">demand from g(x)</small><small className="supply-label">supply Hˢ(r)</small>
        </div>
      </section>
    );
  }

  if (object.area === "Renewal") {
    return (
      <section className="context-graphic renewal-graphic">
        {caption}
        <div className="cohort-flow">
          <div><span>F</span><strong>births</strong></div><i>→</i>
          <div><span>𝓜</span><strong>maturation</strong></div><i>→</i>
          <div><span>B₀</span><strong>local entrants</strong></div><i>+</i>
          <div className="outside-flow"><span>M</span><strong>outside inflow</strong></div><i>→</i>
          <div><span>S</span><strong>stationary scale</strong></div>
        </div>
      </section>
    );
  }

  if (object.area === "Outputs") {
    return (
      <section className="context-graphic policy-graphic">
        {caption}
        <div className="counterfactual-grid">
          <div><span>BASELINE</span><strong>θ, policy₀</strong><small>solve p₀, g₀, S₀, T₀</small></div>
          <i>Δ policy</i>
          <div><span>COUNTERFACTUAL</span><strong>same θ, policy₁</strong><small>resolve p₁, g₁, S₁, T₁</small></div>
        </div>
      </section>
    );
  }

  if (object.area === "Historical context") {
    return (
      <section className="context-graphic historical-graphic">
        {caption}
        <div><span>CENTER</span><i>sorting + prices</i><span>PERIPHERY</span></div>
        <p>Mechanism reference only: its targets, geography, and loss are not current-contract objects.</p>
      </section>
    );
  }

  return null;
}

function relatedObjects(object: MapObject) {
  const ids: string[] = [];
  for (const edge of EDGES) {
    if (edge.from === object.id) ids.push(edge.to);
    if (edge.to === object.id) ids.push(edge.from);
  }
  for (const trace of TRACES) {
    if (trace.nodes.includes(object.id) || trace.packets.includes(object.id)) ids.push(...trace.nodes, ...trace.packets);
  }
  return [...new Set(ids)]
    .filter((id) => id !== object.id && OBJECTS[id])
    .slice(0, 5)
    .map((id) => OBJECTS[id]);
}

export default function Home() {
  const [selectedId, setSelectedId] = useState("target-contract");
  const [traceId, setTraceId] = useState("calibration");
  const [paused, setPaused] = useState(false);
  const [statusFilter, setStatusFilter] = useState<ObjectStatus | null>(null);
  const [copied, setCopied] = useState<string | null>(null);
  const [query, setQuery] = useState("");

  const trace = TRACES.find((item) => item.id === traceId) ?? TRACES[0];
  const selected = OBJECTS[selectedId] ?? OBJECTS["target-contract"];
  const explanation = EXPLANATIONS[selected.id];
  const activeNodes = useMemo(() => new Set(trace.nodes), [trace.nodes]);
  const activeEdges = useMemo(() => new Set(trace.edges), [trace.edges]);

  useEffect(() => {
    function syncFromHash() {
      const hashId = window.location.hash.slice(1);
      if (OBJECTS[hashId]) setSelectedId(hashId);
    }
    window.addEventListener("hashchange", syncFromHash);
    syncFromHash();
    return () => window.removeEventListener("hashchange", syncFromHash);
  }, []);

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
    selectObject(next.nodes[0]);
  }

  function selectObject(nextId: string) {
    setSelectedId(nextId);
    window.history.replaceState(null, "", `${window.location.pathname}${window.location.search}#${nextId}`);
  }

  const prompt = [
    `Investigate the project object: ${selected.title}.`,
    `Economic meaning: ${selected.meaning}`,
    `Plain-English intuition: ${explanation?.intuition ?? selected.meaning}`,
    `Why it matters: ${explanation?.whyItMatters ?? "Verify its role in the equilibrium system."}`,
    `Economic roles: ${(RICH_CONTENT[selected.id]?.roles ?? []).map((role) => role.label).join(", ") || "See the selected object contract."}`,
    `Mechanism steps:\n${(RICH_CONTENT[selected.id]?.steps ?? []).map((item, index) => `${index + 1}. ${item}`).join("\n")}`,
    `Key equation (LaTeX): ${selected.equation}`,
    `Current trace: ${trace.label} — ${trace.question}`,
    `Start from these exact files:\n${selected.sources.map((path) => `- ${path}`).join("\n")}`,
    `Known caveats:\n${selected.caveats.map((item) => `- ${item}`).join("\n")}`,
    `Debug questions:\n${(explanation?.debugQuestions ?? []).map((item) => `- ${item}`).join("\n")}`,
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
          {repositoryLink("CALIBRATION_STATUS.md") && (
            <a href={repositoryLink("CALIBRATION_STATUS.md") ?? undefined} target="_blank" rel="noreferrer">
              Open canonical status ↗
            </a>
          )}
        </div>
      </header>

      <SystemGuide selectObject={selectObject} />

      <ResearchSnapshot selectObject={selectObject} />

      <ObjectFinder query={query} setQuery={setQuery} selectObject={selectObject} />

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

      <TraceRoute trace={trace} selectedId={selectedId} selectObject={selectObject} />

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
                    <button onClick={() => selectObject(node.id)} className={selectedId === node.id ? "selected" : ""}>
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
                        onClick={() => selectObject(packetId)}
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
                      onClick={() => selectObject(packetId)}
                      onKeyDown={(event) => {
                        if (event.key === "Enter" || event.key === " ") selectObject(packetId);
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
                    onClick={() => selectObject(node.id)}
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
                  <button key={packetId} onClick={() => selectObject(packetId)} className={selectedId === packetId ? "selected" : ""}>
                    <i className={statusSlug(packet.status)} /> {packet.title}
                  </button>
                );
              })}
            </div>
          </div>
        </div>

        <Inspector
          key={selected.id}
          object={selected}
          traceLabel={trace.label}
          copied={copied}
          copyText={copyText}
          prompt={prompt}
          selectObject={selectObject}
        />
      </section>

      <footer className="footer-note">
        <span>MAP CONTRACT</span>
        <p>
          Status, target values, weights, bounds, fingerprints, and closure diagnostics are refreshed from canonical repository files. Economic descriptions and graph topology are intentionally curated.
        </p>
        <div className="footer-links">
          <nav aria-label="Canonical project anchors">
            <a href={repositoryLink("CALIBRATION_STATUS.md") ?? undefined} target="_blank" rel="noreferrer">Status</a>
            <a href={repositoryLink("code/model/README.md") ?? undefined} target="_blank" rel="noreferrer">Model index</a>
            <a href={repositoryLink("latex/dynamic_intergenerational_housing_fertility_model.tex") ?? undefined} target="_blank" rel="noreferrer">Paper source</a>
          </nav>
          <code>{statusData.generatedFrom.targetSet} · {statusData.generatedFrom.fingerprint.slice(0, 12)}…</code>
        </div>
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
  selectObject,
}: {
  object: MapObject;
  traceLabel: string;
  copied: string | null;
  copyText: (key: string, text: string) => Promise<void>;
  prompt: string;
  selectObject: (id: string) => void;
}) {
  const [promptOpen, setPromptOpen] = useState(false);
  const explanation = EXPLANATIONS[object.id];
  const rich = RICH_CONTENT[object.id];
  const related = relatedObjects(object);

  return (
    <aside className="inspector" aria-live="polite">
      <div className="inspector-topline">
        <span>{object.kind.toUpperCase()} / {object.area.toUpperCase()}</span>
        <span className={`status-pill ${statusSlug(object.status)}`}>{object.status}</span>
      </div>
      <h2>{object.title}</h2>
      <button
        className="block-link"
        onClick={() => copyText("block-link", `${window.location.origin}${window.location.pathname}#${object.id}`)}
      >
        {copied === "block-link" ? "Block link copied" : "Copy link to this block"} <span>↗</span>
      </button>

      <section className="definition-block">
        <span>DEFINITION</span>
        <p className="meaning">{object.meaning}</p>
      </section>

      {rich && (
        <div className="role-strip" aria-label="Research roles">
          {rich.roles.map((role) => <span key={role.label} className={role.kind}>{role.label}</span>)}
        </div>
      )}

      {explanation && (
        <section className="explanation-card">
          <span>IN PLAIN ENGLISH</span>
          <p>{explanation.intuition}</p>
        </section>
      )}

      <MechanismGraphic object={object} />

      <ContextGraphic object={object} />

      {rich && (
        <section className="mechanism-steps">
          <div><span>HOW THIS BLOCK WORKS</span><small>{rich.steps.length} steps</small></div>
          <ol>
            {rich.steps.map((step, index) => (
              <li key={step}><span>{String(index + 1).padStart(2, "0")}</span><p>{step}</p></li>
            ))}
          </ol>
        </section>
      )}

      <section className="detail-block equation-block">
        <h3>KEY EQUATION</h3>
        <MathExpression tex={object.equation} display />
        {explanation && <p>{explanation.equationReading}</p>}
      </section>

      {explanation && (
        <section className="symbol-glossary">
          <h3>READ THE SYMBOLS</h3>
          <dl>
            {explanation.terms.map((term) => (
              <div key={term.symbol}>
                <dt><MathExpression tex={term.symbol} /></dt>
                <dd>{term.meaning}</dd>
              </div>
            ))}
          </dl>
        </section>
      )}

      {rich?.example && (
        <section className="worked-example">
          <div><span>WORKED EXAMPLE</span><small>illustrative arithmetic</small></div>
          <h3>{rich.example.title}</h3>
          <p>{rich.example.setup}</p>
          <MathExpression tex={rich.example.equation} display />
          <p>{rich.example.result}</p>
        </section>
      )}

      {explanation && (
        <section className="why-block">
          <span>WHY THIS BLOCK MATTERS</span>
          <p>{explanation.whyItMatters}</p>
        </section>
      )}

      {object.id === "target-contract" && <TargetContract />}
      {object.id === "parameters" && <ParameterContract />}
      {object.id === "population-closure" && <ClosureSnapshot />}

      <section className="detail-block sources-block">
        <h3>EXACT SOURCE FILES</h3>
        <ul>
          {object.sources.map((path, index) => {
            const exists = sourceExists(path);
            const remote = repositoryLink(path);
            const canBrowse = Boolean(remote) && !path.startsWith("output/");
            const localPath = `${statusData.generatedFrom.repoRoot}/${path}`;
            return (
              <li key={path}>
                <div className={exists === false ? "source-main missing" : "source-main"}>
                  <i>{exists === false ? "!" : String(index + 1).padStart(2, "0")}</i>
                  <span>
                    <code>{path}</code>
                    <small>{sourceKind(path)} · {exists === false ? "not present in this worktree" : exists ? "verified locally" : "listed source"}</small>
                  </span>
                </div>
                <div className="source-actions">
                  {canBrowse && <a href={remote ?? undefined} target="_blank" rel="noreferrer">View ↗</a>}
                  <button onClick={() => copyText(`source-${index}`, localPath)} aria-label={`Copy absolute path for ${path}`}>
                    {copied === `source-${index}` ? "Copied" : "Copy path"}
                  </button>
                </div>
              </li>
            );
          })}
        </ul>
      </section>

      {explanation && (
        <section className="debug-block">
          <h3>DEBUG THIS OBJECT</h3>
          <ol>
            {explanation.debugQuestions.map((question) => <li key={question}>{question}</li>)}
          </ol>
        </section>
      )}

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

      {related.length > 0 && (
        <section className="related-block">
          <div><span>FOLLOW THE SYSTEM</span><small>related objects</small></div>
          <nav aria-label={`Objects related to ${object.title}`}>
            {related.map((item) => (
              <button key={item.id} onClick={() => selectObject(item.id)}>
                <span>{item.area}</span><strong>{item.title}</strong><i>→</i>
              </button>
            ))}
          </nav>
        </section>
      )}

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
          <thead><tr><th>Moment</th><th>Target</th><th>SE</th><th>Weight</th><th>Basis</th></tr></thead>
          <tbody>
            {statusData.targetRows.map((row) => (
              <tr key={row.name} className={row.status === "empirically held" ? "held-row" : ""}>
                <td>{humanize(row.name)}{row.status === "empirically held" && <small className="held-tag">held</small>}</td>
                <td>{compactNumber(row.target, 5)}</td>
                <td>{compactNumber(row.standardError, 5)}</td>
                <td>{compactNumber(row.weight, 1)}</td>
                <td>{row.weightKind === "measured / declared" ? "measured SE" : "synthetic 5%"}</td>
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
