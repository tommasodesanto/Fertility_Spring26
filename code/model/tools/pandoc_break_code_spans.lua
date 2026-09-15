-- Pandoc Lua filter: allow line breaks inside long code spans and paths so
-- they wrap instead of overflowing the margin. Used for the September 15, 2026
-- structural-review PDFs. Build command:
--   pandoc docs/model/<file>.md --from markdown+tex_math_single_backslash+pipe_tables \
--     --lua-filter=code/model/tools/pandoc_break_code_spans.lua --pdf-engine=xelatex \
--     -V geometry:margin=1in -V fontsize=11pt -V mainfont="Times New Roman" \
--     -V monofont="Menlo" -V colorlinks=true --toc --toc-depth=2 -o output/pdf/<file>.pdf

-- Allow long inline code spans (file paths like docs/model/foo_bar_20260915.md)
-- to break across lines instead of overflowing the page margin.
-- Verified against the three source documents: their inline code spans contain
-- only [A-Za-z0-9] plus / _ . - : , and space, so no generic LaTeX special-
-- character escaping (\ { } $ & % # ^ ~) is needed here.
function Code(el)
  local s = el.text
  s = s:gsub("_", "\\_\\allowbreak{}")
  s = s:gsub("([/.:,-])", "%1\\allowbreak{}")
  return pandoc.RawInline('latex', '\\texttt{' .. s .. '}')
end

-- Some plain-prose Str tokens (not marked up with backticks) are also file
-- paths, e.g. "latex/september_14_presentation.tex:127-145" in the ChatGPT
-- review. Only touch tokens containing "/" (a reliable, low-collision marker
-- for a path or slash-expression like "p90/p50"); everything else is left to
-- pandoc's normal escaping. Verified: no "/"-containing Str token in any of
-- the three source documents contains \ { } $ & % # ^ ~, so only underscore
-- needs escaping here.
function Str(el)
  local s = el.text
  if not s:find("/", 1, true) then return nil end
  s = s:gsub("_", "\\_\\allowbreak{}")
  s = s:gsub("([/.:,-])", "%1\\allowbreak{}")
  return pandoc.RawInline('latex', s)
end
