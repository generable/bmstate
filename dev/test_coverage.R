# Count code coverage of tests
require(covr)
exclusions <- list()
cov <- package_coverage(line_exclusions = exclusions, pre_clean = TRUE)

# HTML report
report(cov)

# upload to codecov.io
# codecov(coverage = cov, token = tok)
