# FormulaComp
Dirty code of some formulas in Biology, and Physics

class Maths:
	def sine(angle):
		sine = math.sin(math.radians(angle))
		return sine

	def cosine(angle):
		cosine = math.cos(math.radians(angle))
		return cosine

	def tangent(angle):
		tangent = math.tan(math.radians(angle))
		return tangent
	def adding(a, b):
		return lambda a, b: a + b

# object = Maths.adding(12, 23)
# print(object(1,2))


class PhysicsKinematics:
	# Tailgating
	def tailgating(iv, reaction_t, gap, acceleration_A, acceleration_B):
		print('TAILGATING')
		print('\n')
		print('GIVEN: \nInitial Velocity: {}\nReaction Time: {}\nGap: {}\nAcceleration of Car A: {}\nAcceleration of Car B: {}\n'.format(iv, reaction_t, gap, acceleration_A, acceleration_B))
		mps = iv/3.6 # converting the velocity in kph to velocity in mps
		print('Velocity in meters per second: ', mps, 'm/s')
		distopA = mps**2/(2*acceleration_A) * (-1)
		print('The stopping distance of object A after deceleration is: ', distopA, 'm')
		distopB = mps**2/(2*acceleration_B) * (-1)
		print('The stopping distance of object B after deceleration is: ', distopB, 'm')
		breactdist = (reaction_t*mps)
		print('Distance covered by object B before deceleration: ', breactdist, 'm')
		print('Total distance travelled by object B: ', distopB + breactdist, 'm')
		safedistance = (distopB + breactdist) - distopA + gap
		print('Safe distance (for two objects not to crash):', safedistance, 'm')
		distcrash = distopA + gap - breactdist # the distance travelled of object B if it were to crash to obj A
		print('Distance of object B till crash: ', distcrash)
		objB_fv_crash = math.sqrt(2*acceleration_B*distcrash + (mps**2))
		print('Final velocity of Object B if crashed: ', objB_fv_crash, 'm/s')

	def freefall(iv, height, t):
		print('\n')
		print('FREE FALL')
		print('\n')
		d = iv*t + (.5*9.8*t**2)
		print('Distance in {} seconds:'.format(t), d ,'m')
		print('\n')

	def projectile(velocity, angle, gravity):
		print('\nPROJECTILE\n')
		# time of flight
		print('--GIVEN--\nInitial VELOCITY: {}\n---Initial ANGLE: {}\n---------GRAVITY: {}\n'.format(velocity,angle,gravity))
		t_overall = 2*velocity * Maths.sine(angle)/gravity
		R = velocity**2*Maths.sine(angle)*2/gravity
		print('Range: ', R, 'm')
		print('Flight Time: ', t_overall, 's')
		
	def centripetal(tanvel, tanvel2, angvel, angvel2, mass, radius, ang_accel, tan_accel, cent_accel, f, T, tacc, an_angle, staticfriction):
		# NOTE THAT THE 'VELOCITY' IS TANGENT VELOCITY, AND 'RADIANPS' IS ANGULAR VELOCITY.
		if T and mass and radius != 0:
			# a 3kg rock attached to a 2-m cord, swings in a horizontal circle so that it makes one revolution in 0.3s. a.) what is the centripetal force on the rock? b.) what is the tangential velocity of the rock? c.) What is the frequency f?
			print('\n--CENTRIPETALSSS--\n\nGIVEN\nTanvel: {}\nTanvel 2: {}\nRadians/s: {}\nRadians/s(final): {}\nMass: {}\nRadius: {}\nAcceleration: {}\nFrequency: {}\nRevolution(T): {}\nTime till acceleration: {}\nTheta: {}\nStatic Friction: {}\n'.format(tanvel, tanvel2, angvel, angvel2, mass, radius, ang_accel, tan_accel, cent_accel, f, T, tacc, an_angle, staticfriction))
			f = 1/.3*60
			circumference = 2*math.pi*radius
			tanvel = circumference/T
			cent_accel = tanvel**2/radius
			cent_accel2 = (4*math.pi**2*radius)/T**2
			ForceC = mass*cent_accel2
			Fc = cent_accel*mass
			centrif = mass*(tanvel**2/radius)
			print('MASS: ', mass)
			print('REVOLUTION TIME: ', T)
			print('RPM or FREQUENCY: ', f)
			print('CIRCUMFERENCE: ', circumference)
			print('TANGENT VELOCITY:', tanvel)
			print('CENTRIPETAL ACCELERATION: ', cent_accel)
			print('CENTRIPERAL ACCELERATION2: ', cent_accel2)
			print('Force Centripetal: ', ForceC)
			print('FORCE CENTRI ma: ', Fc)
			print('FORCE CENTRI mvR:', centrif)
			quit()
		if tanvel and mass and radius != 0:
			# try to find the force that will make the thing move in the said velocity
			print('\n--CENTRIPETALSSS--\n\nGIVEN\nTanvel: {}\nTanvel 2: {}\nRadians/s: {}\nRadians/s(final): {}\nMass: {}\nRadius: {}\nAcceleration: {}\nFrequency: {}\nRevolution(T): {}\nTime till acceleration: {}\nTheta: {}\nStatic Friction: {}\n'.format(tanvel, tanvel2, angvel, angvel2, mass, radius, ang_accel, tan_accel, cent_accel, f, T, tacc, an_angle, staticfriction))
			Fc = mass*tanvel**2/radius
			print('CENTRIPETAL FORCE: ', Fc)
			quit()
		if angvel and angvel2 and tacc != 0:
			#find the angular velocity(omega), the tangential velocity, angular acceleration(alpha), tangential accel.
			print('\n--CENTRIPETALSSS--\n\nGIVEN\nTanvel: {}\nTanvel 2: {}\nRadians/s: {}\nRadians/s(final): {}\nMass: {}\nRadius: {}\nAcceleration: {}\nFrequency: {}\nRevolution(T): {}\nTime till acceleration: {}\nTheta: {}\nStatic Friction: {}\n'.format(tanvel, tanvel2, angvel, angvel2, mass, radius, ang_accel, tan_accel, cent_accel, f, T, tacc, an_angle, staticfriction))
			ang_accel = (radianps2-radianps)/tacc
			tan_accel = (radius(radianps2 - radianps))/tacc

			print('ANGULAR ACCELERATION: ', angular_acceleration)
			print('TANGENT ACCELERATION: ', tangential_acceleration)
			quit()
		if an_angle and radius != 0:
			# getting the raw angular velocity and tangent velocity
			arc_length = radius*an_angle
			angvel = an_angle/arc_length
			tanvel = radius*an_angle
			quit()
		if staticfriction and mass and radius != 0:
			# COMPUTE FOR THE MAXIMUM NUMBER OF REVOLUTIONS, CENTRIPETAL FORCE NORMAL FORCE, F SUB S, AND THE TANGENT VELOCITY OF THE REMOVED OBJECT
			# FORMULAS: Fs = staticfriction*Normalforce
			# Normal force == Centripetal force
			# also Find f
			NormalForce = mass*9.8
			Fs = staticfriction*NormalForce
			# I derived the time of one revolution (T) from the centri accel formula: Ac = 4*math.pi**2*radius/T**2
			# question is, I can't remember where I got the Ac.
			T = math.sqrt((4*math.pi**2*radius)/3.92)
			f = (1/T)*60
			tanvel = (2*math.pi*radius)/T   # I used the formula 4 velocity: (2*pi*R)/T
			Fc = mass*tanvel
			print('FORCE STATIC: ', Fs)
			print('TIME ONE REV (T):',T)
			print('TANGENT VELOCITY:', tanvel)
			print('RPM (f):', f)
			print('CENTRIPETAL FORCE: ', Fc)
			quit()


class Biology:
	def hardyweinberg(p, q, population): # p is the number of dominant allele frequency in the gene pool. q is the recesive allele frequency in the gene pool
	# YOU HAVE TO DETERMINE FIRST THE NUMBER OF P AND Q IN THE GENE POOL.
		print('\nHARDY-WEINBERG PRINCIPLE\n')
		if p and q and population != 0:
			print('--GIVEN--\nDOMINANT: {}\nRECESSIVE: {}\nPOPULATION: {}\n'.format(p,q,population))
			p_pool = p/population 
			q_pool = q/population
			print('P + Q is equals to:', p_pool+q_pool)
			print('Dominants in the Gene Pool:', p_pool)
			print('Recessives in the Gene Pool:', q_pool)
			# what are the odds that you will get a pp, pq, qp, or qq?
			p2 = p_pool*p_pool
			q2 = q_pool*q_pool
			p2q = 2 * p_pool * q_pool
			print('HOMO DOMINANTS LIKELIHOOD (P2):', p2)
			print('HOMO RECESSIVE LIKELIHOOD (Q2):', q2)
			print(' HETEROZYGOUS LIKELIHOOD (2PQ):', p2q)
			hardwein = p2 + p2q + q2
			print('The P2 + 2PQ + Q2:', hardwein)
		if p == 0:
			print('--GIVEN--\nDOMINANT: {}\nRECESSIVE: {}\nPOPULATION: {}\n'.format(p,q,population))
			q_pool = q/population
			p_pool = 1.0 - q_pool
			print('P + Q = ', p_pool + q_pool)
			print('Dominant in Gene Pool: ', p_pool)
			print('Recessive in Gene Pool: ', q_pool)
			# what are the odds that you will get a pp, pq, qp, or qq?
			p2 = p_pool*p_pool
			q2 = q_pool*q_pool
			p2q = 2 * p_pool * q_pool
			print('HOMO DOMINANTS LIKELIHOOD (P2):', p2)
			print('HOMO RECESSIVE LIKELIHOOD (Q2):', q2)
			print(' HETEROZYGOUS LIKELIHOOD (2PQ):', p2q)
			hardwein = p2 + p2q + q2
			print('The P2 + 2PQ + Q2:', hardwein)
		if q == 0:
			print('--GIVEN--\nDOMINANT: {}\nRECESSIVE: {}\nPOPULATION: {}\n'.format(p,q,population))
			p_pool = p/population
			q_pool = 1.0 - p_pool
			print('P + Q = ', p_pool + q_pool)
			print('Dominant in Gene Pool: ', p_pool)
			print('Recessive in Gene Pool: ', q_pool)
			# what are the odds that you will get a pp, pq, qp, or qq?
			p2 = p_pool*p_pool
			q2 = q_pool*q_pool
			p2q = 2 * p_pool * q_pool
			print('HOMO DOMINANTS LIKELIHOOD (P2):', p2)
			print('HOMO RECESSIVE LIKELIHOOD (Q2):', q2)
			print(' HETEROZYGOUS LIKELIHOOD (2PQ):', p2q)
			hardwein = p2 + p2q + q2
			print('The P2 + 2PQ + Q2:', hardwein)
		if population == 0:
			print('--GIVEN--\nDOMINANT: {}\nRECESSIVE: {}\nPOPULATION: {}\n'.format(p,q,population))
			population = p + q
			p_pool = p/population
			q_pool = q/population
			print('P + Q is equals to:', p_pool+q_pool)
			print('Dominants in the Gene Pool:', p_pool)
			print('Recessives in the Gene Pool:', q_pool)
			# what are the odds that you will get a pp, pq, qp, or qq?
			p2 = p_pool*p_pool
			q2 = q_pool*q_pool
			p2q = 2 * p_pool * q_pool
			print('HOMO DOMINANTS LIKELIHOOD (P2):', p2)
			print('HOMO RECESSIVE LIKELIHOOD (Q2):', q2)
			print(' HETEROZYGOUS LIKELIHOOD (2PQ):', p2q)
			hardwein = p2 + p2q + q2
			print('The P2 + 2PQ + Q2:', hardwein)

class Electrics:
	def current(Q,t):
		I = Q/t
		print('CURRENT:', I)

	def electron_mass():
		electronMass = 1.602*10**-7
		print('Mass of ELECTRON:', electronMass)
		

		
# Electrics.current(100,2)
# Electrics.electron_mass()
# PhysicsKinematics.tailgating(iv=100, reaction_t=.45, gap=2.0, acceleration_A=-9.8, acceleration_B=-9.2)

# PhysicsKinematics.freefall(iv=0, height=100, t=3)

# print(dir(random))

# Biology.hardyweinberg(p=500, q=500, population=1000)

# PhysicsKinematics.projectile(velocity=12, angle=30, gravity=9.8)

PhysicsKinematics.centripetal(tanvel=0, tanvel2=0, angvel=0, angvel2=0, mass=0, radius=0, ang_accel=0, tan_accel=0,cent_accel=0,f=0,T=0, tacc=0, an_angle=0, staticfriction=0.4) # FOR THE ACTIVITY

PhysicsKinematics.centripetal(tanvel=0, tanvel2=0, angvel=0, angvel2=0, mass=0, radius=2, ang_accel=0, tan_accel=0, cent_accel=0,f=0, T=0, tacc=0, an_angle=30,staticfriction=0) # FOR PRACTICES

PhysicsKinematics.centripetal(tanvel=0, tanvel2=0, angvel=0, angvel2=0, mass=2, radius=1.5, ang_accel=0, tan_accel=0, cent_accel=0, f=0, T=0, tacc=0, an_angle=0, staticfriction=0)
