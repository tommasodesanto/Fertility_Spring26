export type ObjectExplanation = {
  intuition: string;
  equationReading: string;
  whyItMatters: string;
  terms: { symbol: string; meaning: string }[];
  debugQuestions: string[];
};

export const EXPLANATIONS: Record<string, ObjectExplanation> = {
  "empirical-builders": {
    intuition: "A number becomes usable evidence only after the population, observation unit, timing rule, estimator, and uncertainty are fixed. The builder is the executable definition of that evidence.",
    equationReading: "Moment k is the output of its own estimator applied to raw records under an explicit sample and weighting rule.",
    whyItMatters: "A changed sample or event clock can look like a model improvement even when the economics did not change.",
    terms: [
      { symbol: "\\widehat m_k", meaning: "empirical estimate for moment k" },
      { symbol: "\\mathcal D_k", meaning: "raw data and analysis sample" },
      { symbol: "\\mathcal E_k", meaning: "estimator, fixed effects, weights, and clustering" },
    ],
    debugQuestions: [
      "Does the builder preserve missing outcomes and the intended event window?",
      "Are the sample, fixed effects, weights, and clustered uncertainty recorded next to the estimate?",
    ],
  },
  "provenance-ledger": {
    intuition: "The provenance ledger is the chain of custody for every calibration row. It says not only what the number is, but exactly what economic object the number measures.",
    equationReading: "Each row is a bundle: builder, estimator, sample, fixed effects, clustering, estimate, and uncertainty.",
    whyItMatters: "Without this bundle, a robustness result or stale estimate can silently become a calibration target.",
    terms: [
      { symbol: "\\mathcal P_k", meaning: "provenance record for target k" },
      { symbol: "SE_k", meaning: "sampling uncertainty attached to the exact estimator" },
    ],
    debugQuestions: [
      "Can a reader reproduce the row without consulting a transcript or remembered convention?",
      "Is the calibration source authoritative, or merely the most recent result on disk?",
    ],
  },
  "target-contract": {
    intuition: "The target contract freezes what the model must match before a search begins. It prevents targets or weights from changing opportunistically across candidates.",
    equationReading: "The contract contains twelve empirical moments and their weights; measured standard errors imply inverse-variance weights where available.",
    whyItMatters: "The fingerprint makes two losses comparable only when they were computed against the same signed contract.",
    terms: [
      { symbol: "\\mathcal I_T", meaning: "complete target-and-weight contract" },
      { symbol: "\\widehat m_k", meaning: "target value" },
      { symbol: "w_k", meaning: "objective weight" },
      { symbol: "SE_k", meaning: "standard error or declared uncertainty scale" },
    ],
    debugQuestions: [
      "Does every free parameter have an informative moment or an external restriction?",
      "Does the launcher fingerprint exactly match the profile used by the solver and collector?",
    ],
  },
  "smm-objective": {
    intuition: "SMM asks whether one parameter vector can jointly reproduce the empirical target system after solving the full equilibrium. Every target miss remains visible as its own loss contribution.",
    equationReading: "For each moment, square the model–data gap, multiply by its fixed weight, and add across all rows.",
    whyItMatters: "A low total loss can still conceal one economically central miss, weak identification, or a parameter pushed to a bound.",
    terms: [
      { symbol: "\\theta", meaning: "candidate vector of free structural parameters" },
      { symbol: "m_k(\\theta)", meaning: "model-implied moment after equilibrium is solved" },
      { symbol: "J(\\theta)", meaning: "weighted SMM objective" },
    ],
    debugQuestions: [
      "Which target rows dominate the loss, and what parameter block should they identify?",
      "Was the incumbent evaluated alongside new chains under the identical contract and tolerances?",
    ],
  },
  parameters: {
    intuition: "Parameters translate economic mechanisms into magnitudes: patience shapes saving, child utility shapes fertility, housing preferences shape room demand, and supply parameters shape price incidence.",
    equationReading: "The optimizer may search only inside the declared parameter domain; an owner must also finance the unlevered share of the house value.",
    whyItMatters: "Bounds and external restrictions are part of identification, not implementation details.",
    terms: [
      { symbol: "\\Theta", meaning: "admissible search region" },
      { symbol: "\\phi", meaning: "financed share of owner housing" },
      { symbol: "(1-\\phi)pH", meaning: "down-payment cash requirement" },
    ],
    debugQuestions: [
      "Which parameters are near bounds, and which target variation is meant to identify them?",
      "Are externally fixed objects empirically normalized, theoretically fixed, or merely provisional?",
    ],
  },
  "household-problem": {
    intuition: "At each lifecycle state, the household compares feasible consumption, saving, rental space, ownership, and fertility choices while anticipating survival, income risk, and future housing prices.",
    equationReading: "Current utility plus the survival-weighted continuation value and death-weighted bequest value is maximized over the feasible choice set.",
    whyItMatters: "This is where housing costs can alter fertility through budgets, tenure thresholds, and the value of space rather than through a reduced-form fertility equation.",
    terms: [
      { symbol: "x", meaning: "age, wealth, income, tenure, parity, and child state" },
      { symbol: "s_a", meaning: "survival probability at age a" },
      { symbol: "B(w^e)", meaning: "utility from terminal estate wealth" },
      { symbol: "\\beta", meaning: "discount factor" },
    ],
    debugQuestions: [
      "Are value functions monotone in wealth and continuous around tenure thresholds?",
      "Do Bellman choices and forward-distribution transitions use the same accepted price and fertility probabilities?",
    ],
  },
  "fertility-choice": {
    intuition: "Eligible households compare waiting with attempting the next birth. A logit taste shock smooths that discrete comparison, while fecundity converts an attempt into a realized birth.",
    equationReading: "The attempt probability is the softmax share of the attempt value relative to the wait value, scaled by parity-specific fertility noise.",
    whyItMatters: "Housing affects the hazard only through how it changes the relative lifetime value of another child.",
    terms: [
      { symbol: "V^1(x)", meaning: "value of attempting a birth" },
      { symbol: "V^0(x)", meaning: "value of waiting" },
      { symbol: "\\kappa(n)", meaning: "fertility-choice noise scale at parity n" },
      { symbol: "\\pi^F(x)", meaning: "birth-attempt probability" },
    ],
    debugQuestions: [
      "Are probabilities inside [0,1] and measured at the correct parity and age?",
      "Is a housing-response moment being mistaken for a sequential fertility-hazard target?",
    ],
  },
  "tenure-savings": {
    intuition: "Renters choose a continuous flow of rooms. Owners choose a discrete housing rung, pay the down payment and user cost, and carry liquid and housing wealth into future states and bequests.",
    equationReading: "Ownership is feasible only when cash covers the unfinanced share; terminal utility depends on nonnegative estate wealth.",
    whyItMatters: "The down-payment kink is a central channel from house prices to ownership, saving, and the space available for children.",
    terms: [
      { symbol: "pH", meaning: "market value of the owner housing rung" },
      { symbol: "b'", meaning: "next-period liquid saving" },
      { symbol: "w^e", meaning: "estate wealth at death" },
      { symbol: "\\theta_0,\\theta_1", meaning: "bequest strength and shift" },
    ],
    debugQuestions: [
      "Does the code use (1−phi), not phi, for the down payment?",
      "Are market quantities based on realized tenure rather than conditional renter policies?",
    ],
  },
  kfe: {
    intuition: "The KFE is population accounting. It takes household policy rules and moves mass across wealth, tenure, income, parity, child, age, survival, and exit states until the cross-section repeats.",
    equationReading: "Apply the transition operator to the current distribution, add entrants, and find the fixed distribution that reproduces itself.",
    whyItMatters: "Prices and moments depend on where household mass sits, not only on what one representative policy function looks like.",
    terms: [
      { symbol: "g(x)", meaning: "stationary mass over household states" },
      { symbol: "\\mathcal T_{\\theta,p}", meaning: "policy-induced transition operator" },
      { symbol: "e", meaning: "young entrant mass" },
    ],
    debugQuestions: [
      "Does mass remain nonnegative and add up under survival, births, maturation, and exit?",
      "Which boundary states accumulate suspicious mass: poorest, richest, renter cap, or housing rungs?",
    ],
  },
  "housing-demand": {
    intuition: "Market demand is the number of physical room-equivalents households actually occupy after the tenure choice. It is not utility from housing and not the renter branch evaluated for owners.",
    equationReading: "Integrate renter rooms for realized renters and discrete owner rungs for realized owners over the stationary distribution.",
    whyItMatters: "A quantity-definition error feeds directly into the price-clearing loop and can contaminate every counterfactual.",
    terms: [
      { symbol: "H^D", meaning: "aggregate physical housing demand" },
      { symbol: "1_R,1_O", meaning: "realized renter and owner indicators" },
      { symbol: "h^R(x)", meaning: "renter rooms" },
      { symbol: "H(o'(x))", meaning: "owner rung in room-equivalent units" },
    ],
    debugQuestions: [
      "Are renter and owner quantities expressed in the same room units?",
      "Is the distribution weighting realized choices rather than conditional policy branches?",
    ],
  },
  "price-clearing": {
    intuition: "The price loop guesses an asset price, resolves household behavior and the distribution, measures housing demand, and updates the price implied by the supply curve.",
    equationReading: "Demand relative to the supply normalization determines a target rent; dividing by the user-cost rate converts rent into an asset price.",
    whyItMatters: "A policy’s incidence depends on the endogenous price response, not only on the partial-equilibrium household budget change.",
    terms: [
      { symbol: "H_0", meaning: "housing-supply normalization" },
      { symbol: "\\xi", meaning: "housing-supply elasticity" },
      { symbol: "\\bar r", meaning: "reference rent normalization" },
      { symbol: "u_c", meaning: "user-cost rate converting rent to price" },
    ],
    debugQuestions: [
      "Does the quantity residual clear at the production tolerance?",
      "Is a stationary asset price being confused with a transition-period housing-cost index?",
    ],
  },
  "births-maturation": {
    intuition: "Births do two jobs: they change the current household’s child state and, after maturation, contribute potential young households to the future entrant pool.",
    equationReading: "Aggregate the realized birth hazard over the distribution, then convert the matured child flow into household-equivalent entrants.",
    whyItMatters: "The maturation clock determines when fertility today becomes population renewal tomorrow.",
    terms: [
      { symbol: "F", meaning: "realized birth flow" },
      { symbol: "B_0", meaning: "mature local-born entrant-household flow" },
      { symbol: "m", meaning: "children currently at home" },
    ],
    debugQuestions: [
      "Does parity count children ever born while the child state counts children at home?",
      "Would changing the maturation architecture require re-estimating the target system?",
    ],
  },
  "population-closure": {
    intuition: "The closure reconciles the young households required by the stationary age structure with mature local-born households and an outside inflow. Scale adjusts so those flows balance.",
    equationReading: "At scale S, required entry equals fixed outside inflow plus retained local-born mature entrants; solving the identity gives the stationary scale.",
    whyItMatters: "This separates a household-level fertility response from the population-level effect that depends on migration and retention assumptions.",
    terms: [
      { symbol: "E_0", meaning: "required young entrant households per normalized adult population" },
      { symbol: "B_0", meaning: "mature city-born entrant households" },
      { symbol: "M", meaning: "fixed outside entrant flow" },
      { symbol: "\\rho", meaning: "local-born retention rate" },
      { symbol: "S", meaning: "stationary adult-household scale" },
    ],
    debugQuestions: [
      "Which closure objects are estimated, empirically normalized, externally fixed, or outstanding?",
      "Is the result a stationary comparison or an actual calendar-time transition?",
    ],
  },
  "policy-engine": {
    intuition: "A policy experiment is a new equilibrium contract, not a parameter toggle. It must state the wedge, fiscal treatment, population closure, and which empirical baseline is being held fixed.",
    equationReading: "The counterfactual must simultaneously satisfy household optimality, housing-market clearing, renewal accounting, and any balanced-budget identity.",
    whyItMatters: "Partial-equilibrium or unfunded results can reverse once prices, transfers, and population scale respond.",
    terms: [
      { symbol: "\\tau_H", meaning: "housing or property-tax wedge" },
      { symbol: "T", meaning: "transfer or grant" },
      { symbol: "\\mathcal C", meaning: "declared closure contract" },
    ],
    debugQuestions: [
      "Is every fiscal and population object reconciled with the live production contract?",
      "Does the reported effect survive price clearing and the current empirical hold?",
    ],
  },
  "diagnostics-paper": {
    intuition: "Diagnostics expose the full equilibrium shape before the paper compresses it into tables and claims. Stable figures make successive calibrations comparable.",
    equationReading: "A reproducible artifact maps one solution cache and its contract metadata into a fixed diagnostic set.",
    whyItMatters: "Scalar moments can look acceptable while policies are jagged, boundary mass is pathological, or markets fail locally.",
    terms: [
      { symbol: "\\mathcal A", meaning: "rendered diagnostic artifact" },
      { symbol: "\\mathcal S", meaning: "solution cache" },
      { symbol: "\\mathcal M", meaning: "target, solver, and closure metadata" },
    ],
    debugQuestions: [
      "Do the standard policy, price, quantity, residual, and boundary plots regenerate from one command?",
      "Does paper language distinguish established results from diagnostic or proposed architecture?",
    ],
  },
  "spatial-predecessor": {
    intuition: "The earlier center–periphery model contains explicit location choice and local housing markets. It remains useful for spatial mechanisms but belongs to a different quantitative contract.",
    equationReading: "Each location has its own housing-demand and supply-clearing condition, linked by household sorting across locations.",
    whyItMatters: "Spatial intuition can be borrowed; losses, targets, and normalized quantities cannot be mixed with the retained one-market strand.",
    terms: [
      { symbol: "i\\in\\{C,P\\}", meaning: "center or periphery location" },
      { symbol: "p_i", meaning: "location-specific house price" },
      { symbol: "g_i", meaning: "household distribution in location i" },
    ],
    debugQuestions: [
      "Is the object being used for mechanism comparison or incorrectly treated as a live benchmark?",
      "Are geography, units, and target definitions comparable before any numerical comparison?",
    ],
  },
  "packet-moments": {
    intuition: "The packet keeps a target value attached to its definition and uncertainty as it enters calibration.",
    equationReading: "A moment travels as value, standard error, and provenance—not as a bare scalar.",
    whyItMatters: "This is the smallest unit that protects empirical meaning.",
    terms: [{ symbol: "(\\widehat m_k,SE_k,\\mathcal P_k)", meaning: "complete empirical-moment packet" }],
    debugQuestions: ["Is any component missing or stale?", "Would another estimator produce the same economic object?"],
  },
  "packet-parameters": {
    intuition: "The candidate packet carries one admissible structural parameter vector into the solver.",
    equationReading: "Ten free coordinates move together; external restrictions remain outside the search vector.",
    whyItMatters: "The model must be evaluated as a joint system, not parameter by parameter.",
    terms: [{ symbol: "\\theta", meaning: "one bounded candidate vector" }],
    debugQuestions: ["Is the candidate inside every bound?", "Are transformed and economic units being confused?"],
  },
  "packet-policies": {
    intuition: "Policy arrays are the household’s state-contingent decisions after solving the Bellman problem.",
    equationReading: "Each state maps to saving, housing, tenure, and fertility rules.",
    whyItMatters: "These rules are the bridge from optimization to population dynamics.",
    terms: [{ symbol: "\\sigma_{\\theta,p}(x)", meaning: "state-contingent household policy rule" }],
    debugQuestions: ["Do rules respect feasibility?", "Are the rules evaluated at the accepted equilibrium price?"],
  },
  "packet-distribution": {
    intuition: "This packet says how much household mass occupies every lifecycle state.",
    equationReading: "The distribution is nonnegative and normalized to one in the maintained stationary composition.",
    whyItMatters: "Aggregation weights policies by actual population mass.",
    terms: [{ symbol: "g(x)", meaning: "normalized stationary household distribution" }],
    debugQuestions: ["Does mass conserve?", "Where does mass pile up at numerical boundaries?"],
  },
  "packet-demand": {
    intuition: "The demand packet is the physical room quantity passed from aggregation to market clearing.",
    equationReading: "Integrate realized housing over the household distribution.",
    whyItMatters: "It is the quantity side of the price fixed point.",
    terms: [{ symbol: "H^D", meaning: "aggregate physical housing demand" }],
    debugQuestions: ["Are all quantities in room-equivalent units?", "Does demand use realized tenure?"],
  },
  "packet-prices": {
    intuition: "The price packet returns from market clearing to every household budget and tenure threshold.",
    equationReading: "The solver dampens the gap between the current and supply-implied target price.",
    whyItMatters: "This feedback converts household demand into general-equilibrium incidence.",
    terms: [{ symbol: "p", meaning: "stationary owner-housing asset price" }],
    debugQuestions: ["Is damping converging rather than cycling?", "Does the final quantity residual meet tolerance?"],
  },
  "packet-births": {
    intuition: "Realized births update current family states and seed the future maturation pipeline.",
    equationReading: "Aggregate state-specific realized birth probabilities over household mass.",
    whyItMatters: "This is the bridge between fertility choice and demographic renewal.",
    terms: [{ symbol: "F", meaning: "aggregate realized birth flow" }],
    debugQuestions: ["Is this a flow rather than completed fertility?", "Does timing match the child-state transition?"],
  },
  "packet-entrants": {
    intuition: "Mature local-born household equivalents combine with outside inflow to form the next entrant cohort.",
    equationReading: "The entrant pool is the sum of outside flow and retained local-born mature flow at the stationary scale.",
    whyItMatters: "It closes the intergenerational accounting loop.",
    terms: [{ symbol: "M+\\rho S B_0", meaning: "available young entrant-household flow" }],
    debugQuestions: ["Are local and outside entrant state distributions justified?", "Is retention empirically disciplined?"],
  },
  "packet-policy": {
    intuition: "The policy packet carries a precisely labeled tax, grant, credit, supply, or retention change together with its closure assumptions.",
    equationReading: "A wedge changes the equilibrium mapping and produces a new fixed point.",
    whyItMatters: "A policy without fiscal and closure metadata is not a reproducible counterfactual.",
    terms: [{ symbol: "\\Delta\\omega", meaning: "declared counterfactual policy wedge" }],
    debugQuestions: ["Is the policy funded?", "Which objects are held fixed across baseline and counterfactual?"],
  },
};
