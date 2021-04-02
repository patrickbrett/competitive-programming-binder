# Returns the length of the longest (strictly) increasing sybsequence of v
# Reverse the input for longest decreasing subsequence.
# Complexity: O(NlogN)
def longest_strictly_increasing_subseq(v):
  s = [0 for _ in v]
  k = 0
  for i in range(len(v)):
    # it = lower_bound(s[0], s[k], v[i])
    pass
    # TODO


# Finds a stable matching with the given preferences.
# mpref[i] list male i's preferred matches in order (highest first).
# fpref lists female preferences in the same format.
# Returns a list of each male's match.
# Complexity: O(N^2)
def stable_matching(mpref, fpref):
  n = len(mpref)
  mpair = [None for _ in range(n)]
  fpair = [None for _ in range(n)]
  p = [0 for _ in range(n)]
  forder = [[None for _ in range(n)] for _ in range(n)]

  for i in range(n):
    for j in range(n):
      forder[i][fpref[i][j]] = j
  
  for i in range(n):
    while mpair[i] < 0:
      p[i] += 1
      w = mpref[i][p[i]]
      m = fpair[w]
      if m is None:
        mpair[i] = w
        fpair[w] = i
      elif forder[w][i] < forder[w][m]:
        mpair[m] = -1
        mpair[i] = w
        fpair[w] = i
        i = m
    
  return mpair