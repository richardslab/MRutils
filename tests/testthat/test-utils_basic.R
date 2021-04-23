

test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("allele replacement works as expected", {
	expect_equal(replace_alleles(list("C", "G"), "C=A,G=T"), c("A", "T"))
	expect_equal(replace_alleles(list("C", "G"), "C=C,G=T"), c("C", "T"))
	expect_equal(replace_alleles(list("C", "G"), "C=A,G=G"), c("A", "G"))
	expect_equal(replace_alleles(list("C", "G"), "C=T,G=A"), c("T", "A"))
	expect_equal(replace_alleles(list("C", "T"), "C=A,T=T"), c("A", "T"))
	expect_equal(replace_alleles(list("C", "T"), "C=A,T=C"), c("A", "C"))
	expect_equal(replace_alleles(list("C", "T"), "C=T,T=C"), c("T", "C"))
	expect_equal(replace_alleles(list("A", "G"), "A=A,G=A"), c("A", "A")) # actually should produce a warning
})