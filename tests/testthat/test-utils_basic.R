

test_that("can source utils.R",{
	expect_warning(source("../../scripts/utils.R"))
})
options(warn=2)
test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

test_that("can get_rsid_from_position", {
	expect_warning(expect_equal(get_rsid_from_position("9", 125711603, "A", "C") , "rs10760259"))
	expect_equal(get_rsid_from_position("9", 125711603, "C", "A") , "rs10760259")
})


test_that("can get_rsid_from_position with chr", {
	expect_warning(expect_equal(get_rsid_from_position("chr9", 125711603, "A", "C") , "rs10760259"))
	expect_equal(get_rsid_from_position("chr9", 125711603, "C", "A") , "rs10760259")
})

test_that("get_rsid_from_position fails appropriately", {
	expect_equal(get_rsid_from_position("blarg", 125711603, "C", "A"), as.character(NULL))
	expect_equal(get_rsid_from_position("9", -125711603, "C", "A"), as.character(NULL))
	expect_equal(get_rsid_from_position("0", 125711603, "C", "A"), as.character(NULL))
	expect_equal(get_rsid_from_position("9", 0, "A", "C"), as.character(NULL))
	expect_equal(get_rsid_from_position("9", 125711603, "T", "C"), as.character(NULL))
	expect_equal(get_rsid_from_position("9", 125711603, "5", "1"), as.character(NULL))
})
