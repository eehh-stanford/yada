# Test comb, dual, elem, and elemToIndex in combinadics.R

N <- 5
k <- 3

# comb
expect_equal(
  comb(0,N,k),
  c(2,1,0)
)

expect_equal(
  comb(1,N,k),
  c(3,1,0)
)

expect_equal(
  comb(2,N,k),
  c(3,2,0)
)

expect_equal(
  comb(3,N,k),
  c(3,2,1)
)

expect_equal(
  comb(4,N,k),
  c(4,1,0)
)

expect_equal(
  comb(5,N,k),
  c(4,2,0)
)

expect_equal(
  comb(6,N,k),
  c(4,2,1)
)

expect_equal(
  comb(7,N,k),
  c(4,3,0)
)

expect_equal(
  comb(8,N,k),
  c(4,3,1)
)

expect_equal(
  comb(9,N,k),
  c(4,3,2)
)

# dual
for(m in 0:(choose(N,k)-1)) {
  expect_equal(
    dual(m,N=N,k=k),
    choose(N,k)-m-1
  )

  expect_equal(
    dual(m,nElements=choose(N,k)),
    choose(N,k)-m-1
  )
}

# elem
expect_equal(
  elem(0,N,k),
  c(0,1,2)
)

expect_equal(
  elem(1,N,k),
  c(0,1,3)
)

expect_equal(
  elem(2,N,k),
  c(0,1,4)
)

expect_equal(
  elem(3,N,k),
  c(0,2,3)
)

expect_equal(
  elem(4,N,k),
  c(0,2,4)
)

expect_equal(
  elem(5,N,k),
  c(0,3,4)
)

expect_equal(
  elem(6,N,k),
  c(1,2,3)
)

expect_equal(
  elem(7,N,k),
  c(1,2,4)
)

expect_equal(
  elem(8,N,k),
  c(1,3,4)
)

expect_equal(
  elem(9,N,k),
  c(2,3,4)
)

# elemToIndex
expect_equal(
  elemToIndex(c(0,1,2),N),
  0
)

expect_equal(
  elemToIndex(c(0,1,3),N),
  1
)

expect_equal(
  elemToIndex(c(0,1,4),N),
  2
)

expect_equal(
  elemToIndex(c(0,2,3),N),
  3
)

expect_equal(
  elemToIndex(c(0,2,4),N),
  4
)

expect_equal(
  elemToIndex(c(0,3,4),N),
  5
)


expect_equal(
  elemToIndex(c(1,2,3),N),
  6
)

expect_equal(
  elemToIndex(c(1,2,4),N),
  7
)

expect_equal(
  elemToIndex(c(1,3,4),N),
  8
)

expect_equal(
  elemToIndex(c(2,3,4),N),
  9
)
