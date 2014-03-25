import numpy as np

# Here I am just experimenting with inheritance and composition.
# Inheritance is when a second (child or sub class) class gets the attributes of a first (parent or base class) 
# Composition is when you pass a class instance or module itself through to the second class and just call the first class through the second.


class Test1():

	def __init__(self):
		self.initial = 10

	def multiply(self,x,y):
		self.constant = 45
		self.result = x*y

T1 = Test1()

class Test2(object,Test1):

	def __init__(self,T1):
		self.T1 = T1

	def divide(self,x,y):
		self.constant2 = 55
		self.result2 = x/y

T2 = Test2(T1)

##
print ''

print '-'*20
print 'Initial Class Dictionaries'
print 'T1.__dict__	=',T1.__dict__
print 'T2.T1.__dict__	=',T2.T1.__dict__
print 'T2.__dict__	=',T2.__dict__
print '-'*20
print ''

T2.multiply(2.0,3.0)
print '-'*20
print 'After T2.multiply(2.0,3.0) call'
print 'T1.__dict__	=',T1.__dict__
print 'T2.T1.__dict__	=',T2.T1.__dict__
print 'T2.__dict__	=',T2.__dict__
del T2.__dict__['result']
del T2.__dict__['constant']
print 'Revert Back to Inital Dictionaries'
print '-'*20
print ''

T2.T1.multiply(3.0,4.0)
print '-'*20
print 'After T2.T1.multiply(3.0,4.0) call'
print 'T1.__dict__	=',T1.__dict__
print 'T2.T1.__dict__	=',T2.T1.__dict__
print 'T2.__dict__	=',T2.__dict__
del T1.__dict__['result']
del T1.__dict__['constant']
print 'Revert Back to Initial Dictionaries'
print '-'*20
print ''

T2.divide(10.0,5.0)
print '-'*20
print 'After T2.divide(10.0,5.0) call'
print 'T1.__dict__	=',T1.__dict__
print 'T2.T1.__dict__	=',T2.T1.__dict__
print 'T2.__dict__	=',T2.__dict__
del T2.__dict__['result2']
del T2.__dict__['constant2']
print 'Revert Back to Initial Dictionaries'
print '-'*20

