# Log
## [2026-04-08 16:43:38] ingest | channel metadata sync

- channel_url: `https://www.youtube.com/channel/UC2ggjtuuWvxrHHHiaDH1dlQ/videos`
- video_count: `478`
- output: `raw/youtube/channel_videos.json`
## [2026-04-08 16:44:15] compile | wiki refresh

- video_count: `478`
- topic_pages: `8`
- series_pages: `203`
- cached_transcripts: `0`
## [2026-04-08 16:44:15] lint | wiki health check

- problems_found: `1`
- report: `wiki/lint-report.md`
## [2026-04-08 16:47:58] ingest | transcript sync

- selected_videos: `80`
- fetched_now: `27`
- missing_now: `22`
- errors_now: `31`
- languages: `zh-TW, zh-Hant, zh`
## [2026-04-08 16:48:15] compile | wiki refresh

- video_count: `478`
- topic_pages: `8`
- series_pages: `203`
- cached_transcripts: `27`
## [2026-04-08 16:48:15] query | 什麼是 transformer attention

- results: `6`
- output: `outputs/query-briefs/什麼是-transformer-attention.md`
- wiki_copy: `wiki/queries/什麼是-transformer-attention.md`
## [2026-04-08 16:48:15] lint | wiki health check

- problems_found: `1`
- report: `wiki/lint-report.md`
## [2026-04-08 16:48:24] compile | wiki refresh

- video_count: `478`
- topic_pages: `8`
- series_pages: `203`
- cached_transcripts: `27`
## [2026-04-08 16:48:24] lint | wiki health check

- problems_found: `0`
- report: `wiki/lint-report.md`
## [2026-04-08 16:49:09] compile | wiki refresh

- video_count: `478`
- topic_pages: `8`
- series_pages: `203`
- cached_transcripts: `27`
## [2026-04-08 16:49:09] lint | wiki health check

- problems_found: `0`
- report: `wiki/lint-report.md`
## [2026-04-08 16:57:26] query | 什麼是 transformer attention

- results: `4`
- output: `outputs/query-briefs/什麼是-transformer-attention.md`
- wiki_copy: `wiki/queries/什麼是-transformer-attention.md`
## [2026-04-08 16:57:26] compile | wiki refresh

- video_count: `478`
- topic_pages: `8`
- series_pages: `203`
- cached_transcripts: `27`
## [2026-04-08 16:57:26] lint | wiki health check

- problems_found: `0`
- report: `wiki/lint-report.md`
## [2026-04-08 16:58:30] query | 什麼是 context engineering

- results: `4`
- output: `outputs/query-briefs/什麼是-context-engineering.md`
- wiki_copy: `wiki/queries/什麼是-context-engineering.md`
## [2026-04-08 16:58:31] compile | wiki refresh

- video_count: `478`
- topic_pages: `8`
- series_pages: `203`
- cached_transcripts: `27`
## [2026-04-08 23:44:00] schema | SKILL.md and persona deep overhaul

- action: Rewrote `SKILL.md`, `references/persona.md`, and `wiki/teaching-style.md`
- reason: Deeper transcript analysis to capture Hung-Yi Lee's authentic teaching patterns
- key additions:
  - 10 named rhetorical moves with transcript-backed examples
  - 「你可能會想說」anticipation pattern
  - 「一言以蔽之」one-sentence-punch pattern
  - Scale-and-surprise with concrete numbers
  - Transition marker table with Chinese phrases
  - Vivid metaphor catalog (暗房裡的人, 擲骰子, 文字接龍, 餓狼下坡)
  - Expanded default response shape (10-step)
  - Paper explanation template
  - Debugging workflow with symptom-first approach
  - Enriched topic priorities with AI Agent and LLM architecture details
## [2026-04-08 23:52:00] schema | spirit deep-dive and SKILL.md philosophy integration

- action: Created `references/spirit.md`, updated `SKILL.md` with Teaching Spirit section and structured guardrails
- reason: Deeper transcript analysis focusing on philosophical values, not just rhetorical patterns
- transcripts analyzed in depth: TigfpYPJk1s (GenAI intro), Taj1eHmZyWw (ML fundamentals), dWQVY_h0YXU (evaluation pitfalls), YJoegm7kiUM (LLM learning journey), bJFtcwLSNxI (DeepSeek R1 reasoning), s266BzGNKKc (LLM evaluation issues), 2rcJdFuNbZQ (AI agents)
- key spirit principles added:
  - Intellectual Honesty First — say when things are hard, uncertain, or heuristic
  - Scale Demystification — make numbers tangible (15T tokens = 1500km of A4 paper)
  - Benchmark Skepticism — Goodhart's Law, Parrot experiment, leaderboard contamination
  - Progressive Formalism — name → intuition → formula → code, never formula-first
  - Analogy Lifecycle — introduce, stretch, break, formalize
  - Research as Living Process — papers-as-data-points, geological time scale (上古時代, 寒武紀)
  - Celebrating The Absurd — HuggingFace origin, Microwave GAN, NoClaw
  - Structured guardrails: Honesty, Metric, Style, Analogy categories
## [2026-04-09 00:22:00] schema | fix stiff output for out-of-domain topics

- action: Updated `SKILL.md` and `wiki/query-playbook.md`
- reason: When using the skill to analyze a system card (Claude Mythos), the output lost the teaching voice entirely — no 「你可能會想說…」, no roadmap, no analogies, switching to blog-post/analyst style mid-answer
- root cause: When To Use was scoped too narrowly to ML teaching; no response shape existed for analytical/commentary tasks; query playbook had no fallback for out-of-KB topics
- fixes applied:
  - When To Use: added out-of-domain triggers (report analysis, explicit skill invocation)
  - Tone Persistence rule: once activated, teaching voice must persist to end of response
  - New response shape: Analyze A Report, System Card, Or News (6-step template)
  - Anti-Regression Guardrails: explicit list of patterns to avoid (menu branching, checklists, bolded taglines, analyst tone, Insight blocks)
  - Query Playbook Out-of-KB Fallback: maintain spirit principles even without transcript evidence
## [2026-04-09 00:45:00] schema | add flavor layer — voice rhythm, humor DNA, simplification instinct

- action: Updated `SKILL.md` and `references/spirit.md`
- reason: After structural fix, output still read like a competent analyst with Chinese transition phrases — missing the actual personality (short sentences, oral particles, genuine reactions, humor mechanisms, 「其實就是」demystification)
- key additions:
  - SKILL.md: Voice Rhythm And Flavor section with 5 subsections (Short Sentence Rhythm, Simplification Instinct, Genuine Reactions, Deadpan Absurd, 「其實就是」Demystification) — each with ❌/✅ before/after examples
  - SKILL.md: Report Analysis response shape rewritten to mandate flavor at every step with concrete examples
  - spirit.md: Humor Mechanisms with 5 named humor patterns from transcripts (casual bewilderment, exaggerated precision, mundane comparison, genuine surprise, blunt honesty)
  - spirit.md: Voice As Personality section — short sentence rhythm, self-answering questions, oral particles (喔嘛啊耶欸), 「其實就是」as jargon-buster
## [2026-04-09 21:23:12] graph | knowledge graph build

- nodes: `916`
- edges: `3664`
- communities: `10`
- files_processed: `490`
- output: `wiki/graph`
## [2026-04-09 21:23:12] compile | wiki refresh

- video_count: `478`
- topic_pages: `8`
- series_pages: `203`
- cached_transcripts: `27`
## [2026-04-09 21:30:09] lint | wiki health check

- problems_found: `0`
- report: `wiki/lint-report.md`
## [2026-06-11 00:00:00] schema | signature verbal habits mining pass

- action: Updated `SKILL.md`, `references/persona.md`, `wiki/teaching-style.md`
- reason: A frequency-mining pass over the 27 cached transcripts surfaced high-frequency signature habits the skill had not yet captured; goal is a more authentic Hung-Yi Lee voice
- method: grep-counted candidate phrases across `raw/youtube/transcripts/*.md`, then sampled real context lines to verify usage before documenting
- key findings (occurrences): 比如說 609 (vs 舉例來說 43), 假設 518, 也許 249, 這樣子 230, 而已 160, 等一下 145, 你會發現 135, 神奇 105, 怎麼辦 70, 所謂 66, 就結束了 36, 對不對 27, 莫名其妙 17, 號稱 16, 硬 train 1
- key additions:
  - SKILL.md: Signature Verbal Habits table (frequency-verified, with usage rules)
  - SKILL.md: new core moves — 「怎麼辦呢？」Problem Before Method, 「你會發現」Guided Discovery
  - SKILL.md: 「而已」Deflation Suffix and 「就結束了」Anticlimax Ending flavor subsections
  - SKILL.md: The Sharing Frame 「跟大家分享」 — sharer positioning, not authority
  - SKILL.md: verified opening variants, time-boxing, prerequisite declaration (Phase 0/1); verified sign-offs incl. 「不虛此行」 and deferred-depth bridge (Phase 7)
  - SKILL.md: Default Response Shape expanded to 11 steps with problem-driven pivot; evaluation checklist gained 4 recommended items
  - persona.md: sections 11.5–11.7 (habit frequency table, 怎麼辦 pivot, verified closings)
  - teaching-style.md: marker table expanded from 10 to 40+ entries with the 比如說/舉例來說 14× ratio note
- correction (user feedback): promoted 「硬 train 一發」 from a single transcript variant to a signature catchphrase. The 27 cached transcripts (2024–25 LLM/Agent era) under-represent it, but per the user it is one of the most iconic Hung-Yi Lee phrases. Added a dedicated 「硬 train 一發」 flavor section to SKILL.md (deep-learning-beats-hand-engineering narrative arc + honest inversion「硬 train 一發 train 不起來」), updated both habit tables.
## [2026-06-11 08:30:00] schema | first-person interview protocol

- action: Created `references/interview-protocol.md`
- reason: The user may gain direct access to Hung-Yi Lee; the skill is currently reverse-engineered from output only. An interview can capture what transcripts cannot: direct preference feedback on skill outputs (出戲點標記), negative space (abandoned techniques), lecture-prep decision process, non-lecture registers (office hour / email / meeting), catchphrase provenance, identity-boundary authorization, and his own lecture-evaluation rubric
- structure: 7 parts with concrete questions, time-budget variants (30/60/90 min), interviewing tips (record his rewrites verbatim; keep self-report vs observed-behavior conflicts as data), and a post-interview mapping table from each answer type to the skill file it updates
## [2026-06-11 09:00:00] schema | interview protocol rewritten for strict Q&A format

- action: Rewrote `references/interview-protocol.md`
- reason: User constraint — no A/B comparison materials, no document-marking exercise; the interview is strict one-question-one-answer
- adaptations: preference feedback converted to 引句反應題 (quote read verbatim inside the question, one quote per question); A/B comparison replaced by 反向生成題 (elicit the teacher's own one-sentence explanation, diff against skill output offline — elicitation beats discrimination); questions flattened into a strictly priority-ordered list (核心十題 with ★ must-ask markers + themed 加時題庫 Q11–Q28); marking exercise demoted to optional leave-behind materials appendix
## [2026-06-16 13:10:00] schema | add 10 grounded A/B comparison questions to interview protocol

- action: Added A/B 對比題 section (AB1–AB10) to `references/interview-protocol.md`
- reason: User requested 10 A/B comparison questions despite the one-question-one-answer constraint — reconciled by noting A/B works in Q&A format when both versions are read aloud within a single question (the teacher picks one)
- method: Workflow (wf_567ebf7b-c5a) — 5 parallel miners extracted genuine stylistic forks from SKILL.md/persona.md/spirit.md, the cached transcripts, and the golden/negative examples; a judge agent selected the 10 most diverse + non-leading pairs (each A and B convey the same content, differ on one dimension, neither a strawman). 6 agents, ~338k tokens
- the 10 dimensions: scale tactic (narrative vs clean conversion), mundane-comparison-for-scary-AI, self-invented narrative analogy vs report facts, genuine-reaction interjections density, opening greeting, sentence rhythm, method-introduction order, closing style, ACG vs everyday analogy, benchmark skepticism
- AB1–AB4 marked ★ high-priority: they test self-invented flourishes (no transcript provenance) — the moves most likely to be inauthentic
- each question annotated with which version the skill currently uses + which skill file gets updated per answer; post-processing table and interviewing tips extended accordingly
## [2026-06-17 00:00:00] schema | ingest first-person interview answers (HIGHEST authority)

- action: Updated `SKILL.md`, `references/persona.md`, `references/spirit.md`, `references/examples/report-analysis-golden.md`
- reason: Received Hung-Yi Lee's actual answers to the interview protocol (Q1–Q28 + AB1/AB3 + 小金 questions). First-person feedback outranks all reverse-engineered guesses.
- new top-authority layer: added `SKILL.md` "First-Person Calibration (本人訪談確認)" section that explicitly outranks everything below it; mirrored as spirit.md section 0
- the three rules (Q3) elevated to top of Evaluation Criteria + spirit.md: (1) 脈絡 not 流水帳, (2) 有梗/punchline as load-bearing core, (3) teach the genealogy of a method not the method
- CONFIRMED by him (kept/strengthened): classroom greeting; genuine reactions incl. 你沒有看錯 (Q11); 其實就是…而已 (Q13); original concrete narrative analogies — printer-intern rated 很好的比喻 (Q12, AB3), concrete detail preferred (本人對小金的話「要講一些具體的內容」); problem-before-method natural+deliberate (Q16); less-is-more (Q15)
- CORRECTED (skill was wrong/over-reaching): drop 「熱騰騰」 (Q1, removed from golden example); prefer 莫名其妙的動漫比喻 over plain everyday comparison for scary facts (Q2, rewrote spirit.md humor mechanism 3 + SKILL.md Genuine Reactions); scale comparison must convey significance not just restate (AB1); analogy theory = pre-load shared content, 芙莉蓮魔族 canonical, AVOID 獵人/黑暗大陸 as dated「老人臭」(Q27, rewrote Technique 7); insight > math and dropping the formula is the higher skill (Q6); calibrate to prerequisites AND interest, may decline to explain (Q8); relevance-first / prompting-before-foundations (Q18, rewrote Phase 2 Motivation)
- new first-person guardrails: never negative-evaluate a specific entity (Q4); no sexual/political jokes (Q7); identity line 「受李宏毅教學風格啟發」 authorized (Q5); factual guardrails — channel not monetized, did not author the textbook (Q28)
- 硬 train 一發 (Q9): confirmed his coinage, still used but declining in agent era; added companion insight 「在 agent 時代，想做什麼比會做什麼重要」
- persona.md: new sections — Non-Lecture Registers (office hour roadmap-but-looser, grad-student challenge stance, layperson 文字接龍), Lecture Preparation, Judging Lecture Quality (audience-relative not viral), Facts About The Persona (incl. 小金 context + its north star「用 AI 的角度做最難的實驗但讓大家聽懂」); rewrote Analogy Policy
## [2026-06-17 01:00:00] docs | README portrait + genealogy story in teacher's voice

- action: Updated `README.md` and `README.en.md`; created `assets/` (with `assets/README.md` note)
- portrait: both READMEs reference `./assets/hung-yi-lee.jpg` at the top with the authorized caption「受李宏毅教學風格啟發」. NOTE: image binary could not be written from the pasted screenshot — user must save the photo to `assets/hung-yi-lee.jpg` and commit it (instructions in assets/README.md)
- new section "這個 skill 是怎麼長出來的 / How This Skill Grew Up": tells the build story via the teacher's own Rule 3 (show how the idea was invented) — v1 transcript-only → limit (survivorship bias) → Step 1 Fable 5 frequency mining → Step 2 Fable 5 multi-agent workflow generating the interview protocol → Step 3 real interview, first-person answers outrank everything
- explicitly credits **Fable 5 (claude-fable-5)** for the optimization passes and generating the interview script, per user request
- includes an "assumed vs actually said" table surfacing the interview corrections
## [2026-06-11 07:24:05] lint | wiki health check

- problems_found: `0`
- report: `wiki/lint-report.md`
## [2026-06-17 16:47:51] lint | wiki health check

- problems_found: `0`
- report: `wiki/lint-report.md`
## [2026-06-17 16:52:21] lint | wiki health check

- problems_found: `0`
- report: `wiki/lint-report.md`
## [2026-06-17 16:56:15] lint | wiki health check

- problems_found: `0`
- report: `wiki/lint-report.md`
