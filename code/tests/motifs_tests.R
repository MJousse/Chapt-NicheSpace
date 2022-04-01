# Test
load("code/functions_motifs.R")

A = matrix(c(0,1,0,0,0,1,0,0,0), nrow = 3) # i <- j <- k
B = matrix(c(0,1,1,0,0,1,0,1,0), nrow = 3) # i <- j <-> k -> i
C = matrix(c(0,1,0,1,0,1,0,0,0), nrow = 3) # xxxi -> j -> k
D = matrix(c(0,1,1,1,0,1,1,1,0), nrow = 3) # i <-> j <-> k
E = matrix(c(0,0,1,1,0,1,0,0,0), nrow = 3) # i <- k -> j <- i
G = matrix(c(0,1,1,1,0,1,0,0,0), nrow = 3) # xxxi -> j -> k
H = matrix(c(0,0,1,1,0,1,1,1,0), nrow = 3) # i -> j -> k

test_that("example", {
  expect_equal(which(motif_role(A)$motif_count==1), 5)
  expect_equal(as.integer(which(motif_role(A)$position_count[1,]==1)), 1)
  expect_equal(as.integer(which(motif_role(A)$position_count[2,]==1)), 2)
  expect_equal(as.integer(which(motif_role(A)$position_count[3,]==1)), 3)
  expect_equal(which(motif_role(B)$motif_count==1), 14)
  expect_equal(as.integer(which(motif_role(B)$position_count[1,]==1)), 15)
  expect_equal(as.integer(which(motif_role(B)$position_count[2,]==1)), 14)
  expect_equal(as.integer(which(motif_role(B)$position_count[3,]==1)), 14)
  expect_equal(which(motif_role(C)$motif_count==1), 6)
  expect_equal(as.integer(which(motif_role(C)$position_count[1,]==1)), 16)
  expect_equal(as.integer(which(motif_role(C)$position_count[2,]==1)), 18)
  expect_equal(as.integer(which(motif_role(C)$position_count[3,]==1)), 17)
  expect_equal(which(motif_role(D)$motif_count==1), 16)
  expect_equal(as.integer(which(motif_role(D)$position_count[1,]==1)), 25)
  expect_equal(as.integer(which(motif_role(D)$position_count[2,]==1)), 25)
  expect_equal(as.integer(which(motif_role(D)$position_count[3,]==1)), 25)
  expect_equal(which(motif_role(E)$motif_count==1), 8)
  expect_equal(as.integer(which(motif_role(E)$position_count[1,]==1)), 6)
  expect_equal(as.integer(which(motif_role(E)$position_count[2,]==1)), 4)
  expect_equal(as.integer(which(motif_role(E)$position_count[3,]==1)), 5)
  expect_equal(which(motif_role(G)$motif_count==1), 9)
  expect_equal(as.integer(which(motif_role(G)$position_count[1,]==1)), 12)
  expect_equal(as.integer(which(motif_role(G)$position_count[2,]==1)), 12)
  expect_equal(as.integer(which(motif_role(G)$position_count[3,]==1)), 13)
  expect_equal(which(motif_role(H)$motif_count==1), 15)
  expect_equal(as.integer(which(motif_role(H)$position_count[1,]==1)), 28)
  expect_equal(as.integer(which(motif_role(H)$position_count[2,]==1)), 26)
  expect_equal(as.integer(which(motif_role(H)$position_count[3,]==1)), 27)
})
