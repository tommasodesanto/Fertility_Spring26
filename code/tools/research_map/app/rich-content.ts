export type ResearchRoleKind =
  | "observed"
  | "estimated"
  | "external"
  | "endogenous"
  | "accounting"
  | "diagnostic"
  | "historical"
  | "held";

export type RichObjectContent = {
  roles: { label: string; kind: ResearchRoleKind }[];
  steps: string[];
  example?: {
    title: string;
    setup: string;
    equation: string;
    result: string;
  };
};

export const RICH_CONTENT: Record<string, RichObjectContent> = {
  "empirical-builders": {
    roles: [
      { label: "empirically measured", kind: "observed" },
      { label: "sample-dependent", kind: "diagnostic" },
      { label: "rooms row on hold", kind: "held" },
    ],
    steps: [
      "Define the population, observation window, geography, and unit of analysis before computing a statistic.",
      "Apply the estimator, survey weights, fixed effects, and clustering rule that give the moment its empirical meaning.",
      "Write the estimate, uncertainty, and full provenance record as one indivisible output.",
    ],
    example: {
      title: "Why the estimator is part of the target",
      setup: "The same raw rooms observations can answer different questions with or without household fixed effects.",
      equation: String.raw`\widehat m^{\mathrm{levels}}\neq\widehat m^{\mathrm{HH\,FE}}`,
      result: "A scalar alone is not a target: the sample and estimator must travel with it. This is an illustration, not a new empirical result.",
    },
  },
  "provenance-ledger": {
    roles: [
      { label: "empirical contract", kind: "observed" },
      { label: "author decision", kind: "external" },
      { label: "incomplete upstream metadata", kind: "held" },
    ],
    steps: [
      "Connect each model moment to exactly one authoritative empirical builder.",
      "Record estimator, sample, fixed effects, clustering, uncertainty, and measurement definition.",
      "Block promotion when the scalar exists but the measurement contract is incomplete or disputed.",
    ],
    example: {
      title: "A provenance failure can leave the number unchanged",
      setup: "Two files may both contain 0.80, yet only one may use the approved sample and fixed effects.",
      equation: String.raw`(\widehat m,\mathcal P_a)\neq(\widehat m,\mathcal P_b)`,
      result: "Matching scalars do not establish matching empirical objects; provenance is part of the identity.",
    },
  },
  "target-contract": {
    roles: [
      { label: "empirical discipline", kind: "observed" },
      { label: "fixed during search", kind: "external" },
      { label: "one held row", kind: "held" },
    ],
    steps: [
      "Freeze the target names, values, uncertainty scales, and ordering before any candidate is solved.",
      "Hash the complete contract so launchers and collectors can reject drift or mixed chains.",
      "Check identification: free parameters need informative moments or declared external restrictions.",
    ],
    example: {
      title: "How uncertainty becomes objective weight",
      setup: "Suppose an illustrative target has standard error 0.10.",
      equation: String.raw`SE=0.10\quad\Longrightarrow\quad w=SE^{-2}=100`,
      result: "A one-standard-error miss then contributes one unit to the quadratic objective. This arithmetic is illustrative.",
    },
  },
  "smm-objective": {
    roles: [
      { label: "estimated parameter map", kind: "estimated" },
      { label: "model-generated moments", kind: "endogenous" },
      { label: "numerical search", kind: "diagnostic" },
    ],
    steps: [
      "Transform a candidate vector into economically valid parameter values inside the declared bounds.",
      "Solve the household, distribution, and market fixed point before measuring model moments.",
      "Return every gap and contribution—not only the scalar loss—to the optimizer and collector.",
    ],
    example: {
      title: "One row of the loss",
      setup: "Take an illustrative target 1.90, model moment 1.80, and standard error 0.10.",
      equation: String.raw`\left(\frac{1.80-1.90}{0.10}\right)^2=1`,
      result: "The row contributes one. A reported total loss is interpretable only with all of its row contributions.",
    },
  },
  parameters: {
    roles: [
      { label: "free parameters", kind: "estimated" },
      { label: "external restrictions", kind: "external" },
      { label: "economic primitives", kind: "accounting" },
    ],
    steps: [
      "Separate estimated parameters from empirical normalizations and externally fixed primitives.",
      "Map optimizer coordinates into the model domain using the declared transform and bound.",
      "Inspect whether estimates sit near bounds and which moments identify each parameter block.",
    ],
    example: {
      title: "Financed share versus down payment",
      setup: "If the financed share is 0.80 and a house costs 100, the buyer must bring the remainder.",
      equation: String.raw`(1-\phi)pH=(1-0.80)\times100=20`,
      result: "The threshold is 20, not 80. This example encodes the repository’s discrete-time convention.",
    },
  },
  "household-problem": {
    roles: [
      { label: "dynamic choice", kind: "endogenous" },
      { label: "lifecycle timing", kind: "accounting" },
      { label: "prices taken as given", kind: "external" },
    ],
    steps: [
      "Enter a period with age, assets, income, tenure, parity, and children-at-home states.",
      "Compare feasible fertility, tenure, rooms, consumption, and saving choices at current prices.",
      "Carry the chosen balance sheet and demographic transition into the next age or into a bequest.",
    ],
    example: {
      title: "Why a policy cannot be read from one budget line",
      setup: "A grant relaxes today’s purchase budget, but it also changes saving, future tenure, births, and equilibrium prices.",
      equation: String.raw`\Delta V=\Delta u_0+\beta\,\mathbb E[\Delta V_1]+\Delta V^{\mathrm{price}}`,
      result: "The first-round transfer is only one term in the lifecycle response; the full value problem must be resolved.",
    },
  },
  "fertility-choice": {
    roles: [
      { label: "endogenous hazard", kind: "endogenous" },
      { label: "probabilistic choice", kind: "accounting" },
      { label: "shared-clock architecture", kind: "external" },
    ],
    steps: [
      "Compute the value of waiting and the value of attempting a birth at the same household state.",
      "Translate their value difference into an attempt probability using the fertility taste scale.",
      "Send realized births—not completed fertility directly—into the demographic transition.",
    ],
    example: {
      title: "Value differences become hazards",
      setup: "If attempting is worth 0.20 more than waiting and the logit scale is 0.10, the illustrative attempt probability is high.",
      equation: String.raw`\pi^F=\frac{e^{2}}{1+e^{2}}\approx0.881`,
      result: "The scale controls how sharply households respond to value differences; this is not a calibrated probability.",
    },
  },
  "tenure-savings": {
    roles: [
      { label: "endogenous tenure", kind: "endogenous" },
      { label: "endogenous saving", kind: "endogenous" },
      { label: "estate accounting", kind: "accounting" },
    ],
    steps: [
      "Screen owner choices using the discrete-time down-payment threshold.",
      "Compare renter housing services with owner products using the full budget and continuation value.",
      "Carry liquid assets, housing equity, sale costs, and debt consistently into next-period wealth or bequests.",
    ],
    example: {
      title: "Conditional renter space is not realized housing",
      setup: "A household can have a well-defined renter optimum even when the tenure choice selects ownership.",
      equation: String.raw`h^{\mathrm{realized}}=\mathbf 1_R h^R+\mathbf 1_O H(o')`,
      result: "Aggregation must use the chosen-tenure object, not the renter policy evaluated off path.",
    },
  },
  kfe: {
    roles: [
      { label: "endogenous distribution", kind: "endogenous" },
      { label: "forward accounting", kind: "accounting" },
      { label: "stationary object", kind: "diagnostic" },
    ],
    steps: [
      "Push mass through survival, income, demographic, tenure, asset, and aging transitions.",
      "Add entrant mass in the approved entrant states and remove death mass consistently.",
      "Iterate until the distribution reproduces itself under the same household policies and prices.",
    ],
    example: {
      title: "Stationarity balances outgoing and incoming mass",
      setup: "If 90 of 100 household units remain after transitions, 10 entrant units restore the same total mass.",
      equation: String.raw`100\xrightarrow{\mathcal T}90\quad+\quad e=10\quad=\quad100`,
      result: "This balance is a stationary accounting identity, not a statement that calendar-time population is constant.",
    },
  },
  "housing-demand": {
    roles: [
      { label: "endogenous quantity", kind: "endogenous" },
      { label: "distribution-weighted", kind: "accounting" },
      { label: "physical units", kind: "observed" },
    ],
    steps: [
      "Read each household’s realized housing from its chosen tenure branch.",
      "Weight that quantity by the stationary distribution over all household states.",
      "Scale normalized demand only when the active population closure requires it.",
    ],
    example: {
      title: "Aggregate the chosen branch",
      setup: "With 60 renter households using 3 units and 40 owners using 6 units, illustrative demand is:",
      equation: String.raw`H^D=60\times3+40\times6=420`,
      result: "Using conditional renter policies for all 100 households would measure a different object.",
    },
  },
  "price-clearing": {
    roles: [
      { label: "endogenous price", kind: "endogenous" },
      { label: "market-clearing", kind: "accounting" },
      { label: "supply normalization", kind: "external" },
    ],
    steps: [
      "Guess the asset price and convert it into renter-paid user cost under the current tax convention.",
      "Resolve household behavior and the stationary distribution, then aggregate physical demand.",
      "Update price from the supply curve until the quantity residual satisfies the solver tolerance.",
    ],
    example: {
      title: "A ten-percent demand excess with unit elasticity",
      setup: "If demand is 1.10 times the supply normalization and supply elasticity is one, target rent is 10 percent above its reference.",
      equation: String.raw`H^D/H_0=1.10,\ \xi=1\quad\Rightarrow\quad r^{\mathrm{target}}=1.10\bar r`,
      result: "The asset price follows after dividing by user cost; rent and asset price are not interchangeable labels.",
    },
  },
  "births-maturation": {
    roles: [
      { label: "endogenous birth flow", kind: "endogenous" },
      { label: "cohort timing", kind: "accounting" },
      { label: "architecture provisional", kind: "diagnostic" },
    ],
    steps: [
      "Aggregate realized birth probabilities over the current household distribution.",
      "Advance children through the maintained maturation clock without confusing parity with children at home.",
      "Map mature city-born children into potential entrant households at the model’s entry age.",
    ],
    example: {
      title: "Births and entrants are separated in time",
      setup: "One hundred households with an average birth flow of 0.20 create 20 births now, not 20 entrants now.",
      equation: String.raw`F_t=100\times0.20=20,\qquad B_{0,t+j}=\mathcal M_j(F_t)`,
      result: "The maturation operator determines when and how many of those births enter the future household pool.",
    },
  },
  "population-closure": {
    roles: [
      { label: "stationary scale", kind: "accounting" },
      { label: "fixed outside inflow", kind: "external" },
      { label: "diagnostic closure", kind: "diagnostic" },
    ],
    steps: [
      "Measure required young entrant households and mature local-born entrants in common household units.",
      "Add the approved outside inflow and retention assumption to the renewal identity.",
      "Solve the stationary scale jointly with housing demand and price; do not label the result a transition path.",
    ],
    example: {
      title: "A closure that nests normalized scale",
      setup: "With required entry 0.06, mature local entry 0.05, fixed inflow 0.01, and full retention:",
      equation: String.raw`S=\frac{0.01}{0.06-0.05}=1`,
      result: "The normalization is then an open stationary point. This arithmetic is illustrative, not an empirical migration estimate.",
    },
  },
  "policy-engine": {
    roles: [
      { label: "exogenous policy wedge", kind: "external" },
      { label: "general equilibrium", kind: "endogenous" },
      { label: "counterfactual contract", kind: "diagnostic" },
    ],
    steps: [
      "Define the tax, grant, credit, or supply change and its financing rule relative to a named baseline.",
      "Resolve household choices, the distribution, housing price, population closure, and any fiscal transfer.",
      "Report household, price, population, and fiscal effects separately with the closure assumptions attached.",
    ],
    example: {
      title: "Partial equilibrium is not the counterfactual",
      setup: "A purchase grant raises buyers’ resources directly, but induced demand can also raise the equilibrium price.",
      equation: String.raw`\Delta V^{GE}=\Delta V^{direct}+\frac{\partial V}{\partial p}\Delta p+\frac{\partial V}{\partial T}\Delta T`,
      result: "The sign and incidence require the resolved price and financing response, not only the statutory grant.",
    },
  },
  "diagnostics-paper": {
    roles: [
      { label: "verified readout", kind: "diagnostic" },
      { label: "paper-facing artifact", kind: "external" },
      { label: "promotion provisional", kind: "held" },
    ],
    steps: [
      "Regenerate the stable policy-function, distribution, price, quantity, and boundary-state diagnostic set.",
      "Reconcile every displayed number with the target, parameter, closure, and policy metadata that produced it.",
      "Promote text or figures only when the empirical and model contracts support the same claim.",
    ],
    example: {
      title: "A scalar fit cannot certify a result",
      setup: "Two runs can share one loss number while using different targets, weights, closures, or measurement units.",
      equation: String.raw`J_a=J_b\ \nRightarrow\ \mathcal I_a=\mathcal I_b`,
      result: "Certification needs the complete fit table, parameter table, fingerprint, residuals, and diagnostic packet.",
    },
  },
  "spatial-predecessor": {
    roles: [
      { label: "historical architecture", kind: "historical" },
      { label: "two-location mechanism", kind: "diagnostic" },
      { label: "non-comparable loss", kind: "held" },
    ],
    steps: [
      "Use the archived center–periphery system to understand the project’s original spatial-sorting mechanism.",
      "Keep its geography, housing-unit normalization, target system, and objective attached to every comparison.",
      "Do not move its parameter estimates or scalar losses into the live one-market contract.",
    ],
    example: {
      title: "Why historical losses cannot be ranked",
      setup: "A two-location loss and a one-market loss sum gaps over different objects and normalizations.",
      equation: String.raw`J^{C/P}(\theta)\not\sim J^{1M}(\theta)`,
      result: "The predecessor is mechanism context, not another candidate in the current optimization.",
    },
  },
};
