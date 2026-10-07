# Three-project documentation conventions

Project 1 (arabidopsis_drought), Project 2 (arabidopsis_root_singlecell) and Project 3
use scripts/ for executable workflows and results/ for analysis outputs.
Public-facing README text uses an English descriptive title, states the research
question first, identifies inputs/provenance, provides a project-root run entry,
separates main analysis from robustness checks, lists outputs and states interpretation limits.
Project 3 follows those existing conventions; Chinese reasoning documents support local learning.

Directory differences should reflect data type, not arbitrary redesign: Project 1
retains FASTQ/alignment/counts/QC folders; Projects 2/3 retain single-cell input and
reference folders. No existing Project 1 or Project 2 path is changed by this reorganization.

Core 1 has a clear current version plus exploration and historical archive because
an input error required revision. Core 2 is not represented as completed or prefilled.
A future analysis can follow the same scripts/results/docs division. Do not rewrite
historical documents into apparently always-correct narratives.
