class point:
	R = 0.0
	z = 0.0
	rho = 0.0
	def __init__( self , z=0.0 , R=0.0 , rho=0.0 ):
		self.z = z
		self.R = R
		self.rho = rho

	def __lt__( self , right ):
		if self.z < right.z:
			return True
		return False

	def __le__( self , right ):
		if self.z <= right.z:
			return True
		return False

	def __gt__( self , right ):
		if self.z > right.z:
			return True
		return False

	def __ge__( self , right ):
		if self.z >= right.z:
			return True
		return False
