from turtle import *
import numpy as np

import turtle
  
ws = turtle.Screen()

def func(i, j):
  turtle.goto(i, j)
  print([i, j])
  turtle.write(str(i)+","+str(j))

ws.onclick(func)
ws.mainloop()
