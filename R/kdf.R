kdf <- function(l, m, d)
{S <- 0
for(i in 0:m){S1 <- 0
for(j in i:(l + m)){
  S1 <- S1 + (-1)^(j - i)*choose(l + m - i, j - i) / (2 * d + j + 1) *
    2^(2 * d + j + 1)
}
S <- S + 2 * choose(m, i) * S1 / (2 * d + i)
}
drop(S)
}
