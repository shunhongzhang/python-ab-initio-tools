#!/usr/bin/python
# operate.py

number = "0123456789."
operater ="+-*/%"

def getoperator(exp):
  if '**' in exp:
    op = '**'       # power operation
  elif "+" in exp:
    op ='+'
  elif "-" in exp:
    op ='-'
  elif '*' in exp:
    op ='*'
  elif "/" in exp:
    op ='/'
  elif "%" in exp:
    op ='%'
  else:
    op = ''
  return op

def getoperand(exp):
  op = getoperator(exp)
  if op == "" :      # no valid operator
    x = 0; y = 0 
  n = exp.index(op)  # index of operator
  x1 = exp[:n]       # first operand before the operator
  if op !='**' :
    x2 = exp[n+1:]   # operand when meeting "**"
  else:
    x2 = exp[n+2:]   # second operand after the operator
  x = float(x1); y = float(x2)
  return (x, y)

def result(exp):
  op = getoperator(exp)
  x,y = getoperand(exp)
  if op =='+' :
    ans = x + y
  elif op =='-' :
    ans = x - y
  elif op =='*' :
    ans = x * y
  elif op =='/' :
    ans = x / y
  elif op =='%' :
    ans = x % y
  elif op =='**' :
    ans = x ** y
  return ans

if __name__=='__main__':
   print result ('3+5')

