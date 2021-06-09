# test-getRSID.R

test_that("can get_rsid_from_position", {
  val <- get_rsid_from_position("9", 125711603, "C", "A")
  expect_equal(val, "rs10760259")
  val <- get_rsid_from_position("9", 125711603, "C", "A")
  expect_equal(val, "rs10760259")
  
})


test_that("can get_rsid_from_position with chr", {
  expect_equal(get_rsid_from_position("chr9", 125711603, "C", "A") ,
               "rs10760259")
})


test_that("get_rsid_from_position warns appropriately", {
  expect_warnings(val <-
                    get_rsid_from_position("chr9", 1257116033, "A", "C"), 2)
  expect_equal(val, NULL)
  
  expect_warnings(val <-
                    get_rsid_from_position("9", 125711603, "T", "C"), 2)
  expect_equal(val, NULL)
})


test_that("get_rsid_from_position errors appropriately", {
  expect_error(get_rsid_from_position("blarg", 125711603, "C", "A"))
  expect_error(get_rsid_from_position("9", 0, "A", "C"))
  expect_error(get_rsid_from_position("9", -125711603, "C", "A"))
  expect_error(get_rsid_from_position("9", 125711603, "5", "1"))
  expect_error(get_rsid_from_position("0", 125711603, "C", "A"))
  
})
