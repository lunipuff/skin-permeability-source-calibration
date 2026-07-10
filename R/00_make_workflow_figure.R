source("R/00_config.R")
# R/00_make_workflow_figure.R
# Generate Figure 1: analysis workflow

required_packages <- c(
	"DiagrammeR",
	"DiagrammeRsvg",
	"rsvg"
)

missing_packages <- required_packages[
	!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
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

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

workflow_graph <- DiagrammeR::grViz("
digraph workflow {

	graph [
		layout = dot,
		rankdir = TB,
		bgcolor = white,
		pad = 0.35,
		nodesep = 0.45,
		ranksep = 0.55
	]

	node [
		shape = rectangle,
		style = 'rounded,filled',
		fillcolor = '#F4F6F8',
		color = '#2F3A45',
		fontname = Helvetica,
		fontsize = 16,
		margin = 0.18
	]

	edge [
		color = '#2F3A45',
		arrowsize = 0.7,
		penwidth = 1.4
	]

	raw [
		label = 'Raw literature-derived\\nskin permeability dataset'
	]

	clean [
		label = 'Data cleaning and\\ncompound identifier harmonization'
	]

	descriptor [
		label = 'Descriptor screening and\\nredundancy analysis'
	]

	loco [
		label = 'Leave-one-compound-out\\ncandidate model search'
	]

	selected [
		label = 'Selected interpretable\\nQSPR model'
		fillcolor = '#E8F0FE'
	]

	benchmarks [
		label = 'Benchmark comparison:\\nclassical, linear, RF, RDKit'
	]

	domain [
		label = 'Applicability-domain\\nanalysis'
	]

	sensitivity [
		label = 'Ablation and\\nvalidation-sensitivity analyses'
	]

	outputs [
		label = 'QSPR skin permeability model\\nwith defined applicability domain'
		fillcolor = '#EAF7EA'
	]

	raw -> clean
	clean -> descriptor
	descriptor -> loco
	loco -> selected
	selected -> benchmarks
	selected -> domain
	selected -> sensitivity
	benchmarks -> outputs
	domain -> outputs
	sensitivity -> outputs
}
")

svg_text <- DiagrammeRsvg::export_svg(workflow_graph)

writeLines(svg_text, path_workflow_svg)

rsvg::rsvg_png(
	charToRaw(svg_text),
	file = path_workflow_png,
	width = 2400,
	height = 3000
)

rsvg::rsvg_pdf(
	charToRaw(svg_text),
	file = path_workflow_pdf,
	width = 8,
	height = 10
)

message("Workflow figure created.")