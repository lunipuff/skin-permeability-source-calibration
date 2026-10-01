palette <- function(n = 100) {
	bay_modified <- c(
		"#00496f",
		"#0f85a0",
		"#20b2aa",
		"#f0ead6",
		"#f1d94a",
		"#ed8b00",
		"#dd4124"
	)

	grDevices::colorRampPalette(
		bay_modified
	)(n)
}
