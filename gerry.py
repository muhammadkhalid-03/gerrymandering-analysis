# gerry.py
# Eric A. Autry
# CSC 395 Spring 2024
# 05/09/24

import numpy as np
import matplotlib.pyplot as plt

"""
isContiguous: will check if the district labelled with the given val is
              contiguous in the given plan.

plan: districting plan matrix
val: district label to check
N: the size of the plan (NxN)

returns: boolean indicating contiguous or not
"""
def isContiguous(plan, val, N):
  # Find nodes with the val we want to check.
  x,y = np.where(plan==val)

  # Create a copy plan to store connection info as we search.
  connect = plan.copy()

  # Create the stack and push any node that has the val.
  stack = []
  stack.append((x[0],y[0]))

  # While the stack is not empty.
  while len(stack) > 0:

    # Pop the next vertex.
    (r,c) = stack.pop()

    # Mark it as visited/connected (-1 to be distinct).
    connect[r,c] = -1

    # Check if the neighbors exist.
    #   Check if the neighbors have the same val.
    #     Push to the stack.
    if r > 0:
      if connect[r-1,c] == val:
        stack.append((r-1,c))
    if r < N-1:
      if connect[r+1,c] == val:
        stack.append((r+1,c))
    if c > 0:
      if connect[r,c-1] == val:
        stack.append((r,c-1))
    if c < N-1:
      if connect[r,c+1] == val:
        stack.append((r,c+1))

  # Flag all vertices that were marked -1, and all vertices that have val.
  # Compare those flags to see if all vertices that have val were marked.
  # Return True if we have marked them all, i.e. they were connected.
  if sum(sum((connect==-1)*(plan==val))) == sum(sum(plan==val)):
    return True

  # Not connected.
  return False

"""
onBorder: checks if a given node is on the boundary separating the two 
          districts in the given plan.

row, col: the indices of the currnet node to check
plan: the current districting plan matrix
N: the size of the plan (NxN)

returns: boolean indicating whether this node is on the district border
"""
def onBorder(row, col, plan, N):

  # Get node value.
  node_val = plan[row,col]

  # Check up neighbor.
  # Return True if in other district.
  if row > 0:
    if plan[row-1,col] != node_val:
      return True

  # Check down neighbor.
  if row < N-1:
    if plan[row+1,col] != node_val:
      return True

  # Check left neighbor.
  if col > 0:
    if plan[row,col-1] != node_val:
      return True

  # Check right neighbor.
  if col < N-1:
    if plan[row,col+1] != node_val:
      return True
    
  # If all were the same, return False, not on border.
  return False

"""
countBorder: counts the number of zeros and ones on the border.

border: the list of border nodes, (row,col) pairs
plan: the districting plan matrix

returns: the number of ones on the border, the number of zeros on the border
"""
def countBorder(border, plan):
  # Set init counts.
  numOnes = 0
  numZeros = 0

  # Loop over the border and count.
  for (brow, bcol) in border:
    if plan[brow, bcol] == 0:
      numZeros += 1
    else:
      numOnes += 1

  return numOnes, numZeros

"""
plotDistrict: will generate a plot of the given districting plan, then displays
              and optionally saves the corresponding plot.

plan: the districting plan of 0s and 1s to plot
filename: for saving the image (no image is saved if filename="")
"""
def plotDistrict(plan, filename="/Users/muhammadkhalid/Desktop/CSC395/homework7/pic.png"):
  fig = plt.figure()
  plt.imshow(plan, cmap='viridis', interpolation='nearest')
  if filename != "":
    fig.savefig(filename)
    plt.close(fig)
  else:
    plt.close(fig)
  return
