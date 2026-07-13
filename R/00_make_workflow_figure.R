############################################################
# 00_make_workflow_figure.R
# Generate Figure 1: cheminformatics workflow
############################################################

source("R/00_config.R")

############################################################
# Required packages
############################################################

required_packages <- c(
	"DiagrammeR",
	"DiagrammeRsvg",
	"rsvg"
)

missing_packages <- required_packages[
	!vapply(
		required_packages,
		requireNamespace,
		logical(1),
		quietly = TRUE
	)
]

if (length(missing_packages) > 0) {
	stop(
		"Missing required packages: ",
		paste(missing_packages, collapse = ", "),
		"\nInstall them with:\ninstall.packages(c(",
		paste(sprintf('"%s"', missing_packages), collapse = ", "),
		"))",
		call. = FALSE
	)
}

############################################################
# Output paths
############################################################

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

if (!exists("path_workflow_svg")) {
	path_workflow_svg <- "figures/figure_workflow.svg"
}

if (!exists("path_workflow_png")) {
	path_workflow_png <- "figures/figure_workflow.png"
}

if (!exists("path_workflow_pdf")) {
	path_workflow_pdf <- "figures/figure_workflow.pdf"
}

############################################################
# Figure settings
############################################################

node_font <- "Helvetica"

############################################################
# Build workflow graph
############################################################

workflow_graph <- DiagrammeR::grViz(
	paste0(
"
digraph workflow {

	graph [
		layout = dot,
		rankdir = LR,
		bgcolor = white,
		pad = 0.25,
		nodesep = 0.42,
		ranksep = 0.55,
		splines = ortho
	]

	node [
		shape = rectangle,
		style = 'rounded,filled',
		color = '#2F3A45',
		penwidth = 1.2,
		fontname = '", node_font, "',
		fontsize = 13,
		margin = 0.15
	]

	edge [
		color = '#2F3A45',
		arrowsize = 0.65,
		penwidth = 1.2,
		fontname = '", node_font, "',
		fontsize = 10
	]

	########################################################
	# Phase labels
	########################################################

	curate_label [
		label = '1. CURATE',
		shape = plaintext,
		fontname = '", node_font, "',
		fontsize = 18,
		fontcolor = '#2F3A45'
	]

	model_label [
		label = '2. MODEL',
		shape = plaintext,
		fontname = '", node_font, "',
		fontsize = 18,
		fontcolor = '#2F3A45'
	]

	audit_label [
		label = '3. AUDIT',
		shape = plaintext,
		fontname = '", node_font, "',
		fontsize = 18,
		fontcolor = '#2F3A45'
	]

	output_label [
		label = '4. OUTPUT',
		shape = plaintext,
		fontname = '", node_font, "',
		fontsize = 18,
		fontcolor = '#2F3A45'
	]

	########################################################
	# CURATE
	########################################################

	data [
		label = 'Curated human skin\\npermeability data\\nKp endpoint + descriptors',
		fillcolor = '#EAF1F8'
	]

	clean [
		label = 'Curation and\\nharmonization\\nCAS identifiers, numeric fields,\\nendpoint completeness',
		fillcolor = '#EAF1F8'
	]

	attrition [
		label = 'Modeling dataset\\n479 raw observations\\n321 cleaned observations\\n140 unique compounds',
		fillcolor = '#EAF1F8',
		penwidth = 1.5
	]

	########################################################
	# MODEL
	########################################################

	descriptor [
		label = 'Descriptor audit\\nredundancy screening\\nGSE-informed solubility check',
		fillcolor = '#ECEBFF'
	]

	loco [
		label = 'LOCO-CV formula search\\ncompound-level validation\\nquadratic, log, and\\ninteraction terms',
		fillcolor = '#ECEBFF'
	]

	selected [
		label = 'Selected interpretable\\nQSPR formula',
		fillcolor = '#DAD7FF',
		penwidth = 1.8
	]

	########################################################
	# AUDIT
	########################################################

	benchmark [
		label = 'Benchmark ladder\\nnull → Potts-Guy → linear\\n→ QSPR → RF → RDKit RF',
		fillcolor = '#FFF3D9'
	]

	domain [
		label = 'Applicability domain\\ndescriptor-space coverage\\nerror by domain',
		fillcolor = '#FFF3D9'
	]

	ablation [
		label = 'Interpretability analyses\\nablation + partial effects\\ninteraction effects',
		fillcolor = '#FFF3D9'
	]

	validation [
		label = 'Validation-design sensitivity\\n5-fold row-wise CV\\nleave-one-reference-out CV',
		fillcolor = '#FFF3D9'
	]

	diagnostics [
		label = 'Supplementary diagnostics\\ncoefficient stability\\nresiduals + uncertainty',
		fillcolor = '#FFF3D9'
	]

	########################################################
	# OUTPUT
	########################################################

	output [
		label = 'Benchmarked interpretable\\nskin permeability QSPR model\\nwith defined applicability domain',
		fillcolor = '#E2F5E8',
		penwidth = 1.8
	]

	########################################################
	# Rank structure
	########################################################

	{ rank = same; curate_label; data; clean; attrition }
	{ rank = same; model_label; descriptor; loco; selected }
	{ rank = same; audit_label; benchmark; domain; ablation; validation; diagnostics }
	{ rank = same; output_label; output }

	########################################################
	# Invisible phase alignment
	########################################################

curate_label [
	label = '1. CURATE',
	shape = rectangle,
	style = 'rounded,filled',
	fillcolor = '#DCEAF7',
	color = '#2F3A45',
	fontcolor = '#2F3A45',
	fontname = 'Helvetica-Bold',
	fontsize = 18,
	margin = 0.16
]

model_label [
	label = '2. MODEL',
	shape = rectangle,
	style = 'rounded,filled',
	fillcolor = '#E5E2FA',
	color = '#2F3A45',
	fontcolor = '#2F3A45',
	fontname = 'Helvetica-Bold',
	fontsize = 18,
	margin = 0.16
]

audit_label [
	label = '3. AUDIT',
	shape = rectangle,
	style = 'rounded,filled',
	fillcolor = '#FFF0CC',
	color = '#2F3A45',
	fontcolor = '#2F3A45',
	fontname = 'Helvetica-Bold',
	fontsize = 18,
	margin = 0.16
]

output_label [
	label = '4. OUTPUT',
	shape = rectangle,
	style = 'rounded,filled',
	fillcolor = '#DDF2E3',
	color = '#2F3A45',
	fontcolor = '#2F3A45',
	fontname = 'Helvetica-Bold',
	fontsize = 18,
	margin = 0.16
]

	########################################################
	# Main workflow
	########################################################

	data -> clean
	clean -> attrition
	attrition -> descriptor
	descriptor -> loco
	loco -> selected

	selected -> benchmark
	selected -> domain
	selected -> ablation
	selected -> validation
	selected -> diagnostics

	benchmark -> output
	domain -> output
	ablation -> output
	validation -> output
	diagnostics -> output
}
"
	)
)

############################################################
# Export figure
############################################################

svg_text <- DiagrammeRsvg::export_svg(workflow_graph)

writeLines(
	svg_text,
	path_workflow_svg
)

rsvg::rsvg_png(
	charToRaw(svg_text),
	file = path_workflow_png,
	width = 4200,
	height = 2300
)

rsvg::rsvg_pdf(
	charToRaw(svg_text),
	file = path_workflow_pdf,
	width = 12,
	height = 6.5
)

############################################################
# Console message
############################################################

message("Workflow figure created:")
message("  SVG: ", path_workflow_svg)
message("  PNG: ", path_workflow_png)
message("  PDF: ", path_workflow_pdf)