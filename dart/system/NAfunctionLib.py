#!/usr/bin/env python2.7

import re, math
from numpy import *
from Constants import *

def length(u):
    	
	"""Calculates the length of u."""
    	
	return math.sqrt(numpy.dot(u, u))

def normalize(u):
    
    	"""Returns the normalized vector along u."""
    	
	return u/math.sqrt(numpy.dot(u, u))

def cross(u, v):
    	
	"""Cross product of u and v:Cross[u,v] = {-u3 v2 + u2 v3, u3 v1 - u1 v3, -u2 v1 + u1 v2}"""
    	
	return array([ u[1]*v[2] - u[2]*v[1],
                       u[2]*v[0] - u[0]*v[2],
                       u[0]*v[1] - u[1]*v[0] ], float)

def rmatrix(alpha, beta, gamma):
    
	"""Return a rotation matrix based on the Euler angles alpha,beta, and gamma in radians."""

	cosA = math.cos(alpha)
	cosB = math.cos(beta)
	cosG = math.cos(gamma)

	sinA = math.sin(alpha)
	sinB = math.sin(beta)
	sinG = math.sin(gamma)

	R = array(
            [[cosB*cosG, cosG*sinA*sinB-cosA*sinG, cosA*cosG*sinB+sinA*sinG],
             [cosB*sinG, cosA*cosG+sinA*sinB*sinG, cosA*sinB*sinG-cosG*sinA ],
             [-sinB,     cosB*sinA,                cosA*cosB ]], float)

	assert numpy.allclose(numpy.linalg.determinant(R), 1.0)
	
	return R

def rmatrixu(u, theta):
	
	"""Return a rotation matrix caused by a right hand rotation of theta radians around vector u."""
	
	if numpy.allclose(theta, 0.0) or numpy.allclose(numpy.dot(u,u), 0.0):
            return numpy.identity(3, float)

	x, y, z = normalize(u)
	sa = math.sin(theta)
	ca = math.cos(theta)

	R = array(
            [[1.0+(1.0-ca)*(x*x-1.0), -z*sa+(1.0-ca)*x*y,     y*sa+(1.0-ca)*x*z],
             [z*sa+(1.0-ca)*x*y,      1.0+(1.0-ca)*(y*y-1.0), -x*sa+(1.0-ca)*y*z],
             [-y*sa+(1.0-ca)*x*z,     x*sa+(1.0-ca)*y*z,      1.0+(1.0-ca)*(z*z-1.0)]], float)

	try:
            assert numpy.allclose(numpy.linalg.determinant(R), 1.0)
	except AssertionError:
            print "rmatrixu(%s, %f) determinant(R)=%f" % (
        	u, theta, numpy.linalg.determinant(R))
            raise

	return R

def dmatrix(alpha, beta, gamma):
	
	"""Returns the displacement matrix based on rotation about Euler angles alpha, beta, and gamma."""
	
	return rmatrix(alpha, beta, gamma) - numpy.identity(3, float)

def dmatrixu(u, theta):
	
	"""Return a displacement matrix caused by a right hand rotation of theta radians around vector u."""
	
	return rmatrixu(u, theta) - numpy.identity(3, float)

def rmatrixz(vec):
	
	"""Return a rotation matrix which transforms the coordinate system such that the vector vec is aligned along the z axis."""
	u, v, w = normalize(vec)

	d = math.sqrt(u*u + v*v)

	if d != 0.0:
            Rxz = array([ [  u/d, v/d,  0.0 ],
                       	[ -v/d, u/d,  0.0 ],
                       	[  0.0, 0.0,  1.0 ] ], float)
	else:
            Rxz = numpy.identity(3, float)

	Rxz2z = array([ [   w, 0.0,    -d],
                        [ 0.0, 1.0,   0.0],
                        [   d, 0.0,     w] ], float)

	R = matrixmultiply(Rxz2z, Rxz)

	try:
            assert numpy.allclose(numpy.linalg.determinant(R), 1.0)
	except AssertionError:
            print "rmatrixz(%s) determinant(R)=%f" % (vec, numpy.linalg.determinant(R))
            raise

	return R

def calc_distance(a1, a2):
    
    	"""Returns the distance between two argument atoms."""
    	
	if a1 == None or a2 == None:
        	return None
    	
	return length(a1.position - a2.position)

def calc_angle(a1, a2, a3):
	
	"""Return the angle between the three argument atoms."""
	
	if a1 == None or a2 == None or a3 == None:
        	return None
	
	a21 = a1.position - a2.position
	a21 = a21 / (length(a21))

	a23 = a3.position - a2.position
	a23 = a23 / (length(a23))

	return math.acos(numpy.dot(a21, a23))

def calc_torsion_angle(a1, a2, a3, a4):
    	
	"""Calculates the torsion angle between the four argument atoms."""
    	
	if a1 == None or a2 == None or a3 == None or a4 == None:
        	return None

	a12 = a2.position - a1.position
	a23 = a3.position - a2.position
	a34 = a4.position - a3.position

	n12 = cross(a12, a23)
	n34 = cross(a23, a34)

	n12 = n12 / length(n12)
	n34 = n34 / length(n34)

	cross_n12_n34  = cross(n12, n34)
	direction      = cross_n12_n34 * a23
	scalar_product = numpy.dot(n12, n34)

	if scalar_product > 1.0:
            	scalar_product = 1.0
	if scalar_product < -1.0:
        	scalar_product = -1.0

	angle = math.acos(scalar_product)
	
	return angle

def UnitvecToDegree(gltilt,glroll,glangle):
	
	"""Calculate orientation of global angle in Euler angle space"""
	
	ftilt = cmp(gltilt,0)
	froll = cmp(glroll,0)
	
	orient1 = math.degrees(math.acos(pow((gltilt/glangle),2)*ftilt))
	orient2 = math.degrees(math.acos(pow((glroll/glangle),2)*froll))
	
	if ftilt == 1 and froll == 1:
		return orient1
	elif ftilt == -1 and froll == 1:
		return orient1
	elif ftilt == -1 and froll == -1:
		return orient2+90
	elif ftilt == 1 and froll == 0:
		return 0.0		
	else:
		return 360-orient1
	
def DegreeToUnitvec(angle):
	
	"""Extract the contribution of global roll and tilt to the global angle by solving the unit cirkel equation"""

	unitcirkel = []

	unitcirkel.append(math.cos(math.radians(angle)))
	unitcirkel.append(math.sin(math.radians(angle)))

	return unitcirkel

def AngleVector(angle,rcont,tcont):

	"""Construct vector of global roll and tilt from global angle and contribution of roll and tilt to this angle"""

	vector = []

	vector.append((sqrt((pow(angle,2))*abs(rcont)))*cmp(rcont,0))
	vector.append((sqrt((pow(angle,2))*abs(tcont)))*cmp(tcont,0))

	return vector

def AccTwist(refbp,twist):			               

	"""Caclculates the accumulated twist relative to the reference base-pair.
	   The calculation procedure is different depending if the reference plane
	   lies on a base-pair or between two base-pairs
	   * If between base-pairs the basepair at "basenumber" of the reference
	     base-pair is duplicated and inserted at the first position next to the
	     reference base-pair. Each identical base-pair twist values is recalculated
	     accoring to its postion relative to the "fraction" of the reference
	     basepair. The result is one extra value of twist.
	   * If the reference plane coincides with a base-pair no duplication is 
	     nesacary for calculation len(ctwist) == number of base-pairs"""

	base = int(refbp)					
	frup = refbp-base
	frdw = 1-frup

	if frup == 0:						#If reference plane is one a base-pairs
		v1up = twist[0:base]
		v1dw = twist[base:len(twist)]

		v1dw[0] = 0.0

		v1up.reverse()
		v1arrup = array(v1up)
		accup = add.accumulate(v1arrup)

		v1arrdw = array(v1dw)
		accdw = subtract.accumulate(v1arrdw)

		v1up = accup.tolist()
		v1dw = accdw.tolist()
		v1up.reverse()

		return (v1up+v1dw)

	else:							#If reference frame lies between two base-pairs
		v1up = twist[0:base]
		v1dw = twist[base:len(twist)]

		v1up[base-1] = v1up[base-1]*frup
		v1dw[0] = v1dw[0]*frdw

		v1dw.insert(0,0.0)

		v1up.reverse()
		v1arrup = array(v1up)
		accup = add.accumulate(v1arrup)

		v1arrdw = array(v1dw)
		accdw = subtract.accumulate(v1arrdw)

		v1up = accup.tolist()
		v1dw = accdw.tolist()
		v1up.reverse()

		return (v1up+v1dw[0:len(v1dw)])

def TwistCorrect(acctwist,glangles):
		
	"""TwistCorrect corrects the global direction of the bend angle at each base-pair step for
	   the twist angle by adding/substracting values for accumulated twist from the reference 
	   base-pair"""

	corrected = []

	for angle in range(len(glangles)):
		corrected.append(glangles[angle]+acctwist[angle])

	return corrected

def CalculateDistance(vector1,vector2):
	
	"""Calculate the distance between two vectors"""
	
	d= float(0)
	
	if len(vector1) == len(vector2):
		for i in range(len(vector1)):
			d=d+(vector1[i]-vector2[i])**2
		d=math.sqrt(d)
 		return d
	else:
		return None	

def Angle(v1, v2):				
	
	"""Calculate angle of vector from two vectors"""
	
	a = array(v1)
	b = array(v2)

	return sqrt((power(a,2))+(power(b,2)))
			
class ConvertSeq:
	
	"""Nucleic acid sequence conversion tools. Class first converts the input to a list of one letter
	   nucleic acid bases for the template and complementary strand. Input is converted to capital
	   letters. For one letter base-sequences a string is also supported as input. The converted input 
	   can supsequently be exported in a differend format. The options are: bases, base-pairs or 
	   base-pair steps all in either one ore three letter code. Mispaired bases will not be converted 
	   to paired forms."""
	
	def __init__(self,seq1=None,seq2=None):
		
		self.type = ''	
		self.input = []
		
		self._TypeCheck(seq1,seq2)
		
	def _MakeCapital(self,seq):
		
		"""Capitalise letters of sequence input, string or list"""
		
		capital = []
		
		if type(seq) == type(''):
			for n in seq:
				capital.append(n.upper())
		elif type(seq) == type([]):
			for seq in seq:
				capital.append(seq.upper())
		else:
			pass		
		
		return capital
	
	def _TypeCheck(self,seq1,seq2):
		
		"""Check if supplied sequence is base, basepair or basepair step"""
		
		sequence = []
		
		if not seq1 == None:
			sequence.append(self._MakeCapital(seq1))
		if not seq2 == None:
			sequence.append(self._MakeCapital(seq2))
				
		for na in sequence[0]:
			tmp = []
			for n in na:
				tmp.append(n)	
			if '-' in tmp or '/' in tmp:
				if len(tmp) == 3 or len(tmp) == 7:
					self.type = 'pair'
				if len(tmp) == 5 or len(tmp) == 13:
					self.type = 'step'
			elif len(tmp) == 1 or len(tmp) == 3:			
				self.type = 'base'
		
		self._ToBaseOne(sequence)
		
	def _ToBaseOne(self,sequence):
		
		"""Convert everything to two strands of one-letter code bases"""
		
		template = []
		complement = []
		
		if self.type == 'base':			#Convert three to one, make complementary strand
			for resid in sequence[0]:
				if resid in BASELIST3_T:
					template.append(BASELIST1_T[BASELIST3_T.index(resid)])
				elif resid in BASELIST1_T:
					template.append(resid)
			if sequence[1]:
				for resid in sequence[1]:
					if resid in BASELIST3_T:
						complement.append(BASELIST1_T[BASELIST3_T.index(resid)])
					elif resid in BASELIST1_T:
						complement.append(resid)
			else:
				for resid in template:
					complement.append(BASELIST1_C[BASELIST1_T.index(resid)])		
						
		elif self.type == 'pair':		#Convert base-pair three to one. Conserve template and complement
			splitter = re.compile('-')
			single_t = []
			single_c = []
			for resid in sequence[0]:
				single_t.append(splitter.split(resid)[0])
				single_c.append(splitter.split(resid)[1])
			if len(single_t[0]) == 1:
				template = single_t
				complement = single_c
			elif len(single_t[0]) == 3:
				for resid in single_t:
					if resid in BASELIST3_T:
						template.append(BASELIST1_T[BASELIST3_T.index(resid)])
				for resid in single_c:
					if resid in BASELIST3_T:
						complement.append(BASELIST1_T[BASELIST3_T.index(resid)])			
		
		elif self.type == 'step':
			splitter = re.compile('/')
			single_t = []
			single_c = []
			for resid in sequence[0]:
				single_t.append(splitter.split(resid)[0])
				single_c.append(splitter.split(resid)[1])
			
			if len(single_t[0]) == 2:
				count = len(single_t)
				step = 0
				while step < count:
					if step == count-1:
						template.append(single_t[step][0])
						template.append(single_t[step][1])
						complement.append(single_c[step][1])
						complement.append(single_c[step][0])
					else:	
						template.append(single_t[step][0])
						complement.append(single_c[step][1])
					step = step+1
		else:
			template = None
			complement = None		
			
		self.input.append(template)
		self.input.append(complement)
					
	def _ToBaseThree(self,inseq):
		
		"""Convert nucleic acid one letter to three letter code. Input base"""
		
		template = []
		complement = []
		
		for resid in inseq[0]:
			if resid in BASELIST1_T:
				template.append(BASELIST3_T[BASELIST1_T.index(resid)])
		for resid in inseq[1]:
			if resid in BASELIST1_T:
				complement.append(BASELIST3_T[BASELIST1_T.index(resid)])
	
		return template,complement			
	
	def _ToBasepair(self,inseq):
		
		"""Convert nucleic acid base to base-pair. Input list of bases"""

		convert = []

		for na in range(len(inseq[0])):
			pair = []
			pair.append(inseq[0][na])
			pair.append('-')
			pair.append(inseq[1][na])
			pair = ''.join(pair)	
			convert.append(pair)
		
		return convert
		
	def _ToBasepairstep(self,inseq):	
		
		"""Convert nucleic acid base to basepairstep. Input list of bases"""
	
		convert = []
	
		for na in range(len(inseq[0])):
			try:	
				pair = []
				pair.append(inseq[0][na])
				pair.append(inseq[0][na+1])
				pair.append('/')
				pair.append(inseq[1][na+1])
				pair.append(inseq[1][na])
				if len(pair) > 0:
					pair = ''.join(pair)	
					convert.append(pair)
				else:
					convert = None	
			except:
				pass		
		
		return convert
		
	def Export(self,exptype):

		"""Export control module, export options are base1,base3,pair1,pair3,step1,step3"""
		
		if exptype == 'base1':
			return self.input
		elif exptype == 'base3':
			return self._ToBaseThree(self.input)
		elif exptype == 'pair1':
			return self._ToBasepair(self.input)
		elif exptype == 'pair3':
			return self._ToBasepair(self._ToBaseThree(self.input))
		elif exptype == 'step1':
			return self._ToBasepairstep(self.input)
		elif exptype == 'step3':
			return self._ToBasepairstep(self._ToBaseThree(self.input))
		else:
			return None						
		
