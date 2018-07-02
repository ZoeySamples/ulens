import numpy as np
import math
import MulensModel as mm

def get_coefficients(calc, zeta, z1, z2, m, dm):

	zeta_conj = np.conj(zeta)

	if (calc == 'geo_cent'):
		# Specific form of the derived for the geometric center frame

		coeff5 = (z1**2 - zeta_conj**2)

		coeff4 = (-(z1*(2*dm + z1*zeta)) - 2*m*zeta_conj + zeta*zeta_conj**2)

		coeff3 = (-2*z1**4 + 4*(dm*z1 + m*zeta)*zeta_conj + 2*z1**2*zeta_conj**2)

		coeff2 = (4*dm*z1*(m + z1**2) + 2*(2*m**2 + z1**4)*zeta - 4*dm*z1*zeta*zeta_conj - 2*z1**2*zeta*zeta_conj**2)

		coeff1 = z1*(-4*dm**2*z1 - 4*m**2*z1 + z1**5 - 8*dm*m*zeta - 4*z1*(dm*z1 + m*zeta)*zeta_conj - z1**3*zeta_conj**2)

		coeff0 = z1**2*(4*dm*m*z1 - 2*dm*z1**3 + 4*dm**2*zeta - z1**4*zeta + 2*z1*(m*z1 + 2*dm*zeta)*zeta_conj + z1**2*zeta*zeta_conj**2)

	elif (calc == 'plan'):
		# Specific form of the derived for the planet frame

		csum = mm.Utils.complex_fsum

		"""
		coeff5 = - (zeta_conj*(-z2 + zeta_conj))

		coeff4 = (dm*z2 + m*(z2 - 2*zeta_conj) - (2*z2 + zeta)*(z2 - zeta_conj)*zeta_conj)

		coeff3 = (-(dm*z2*(z2 + 2*zeta_conj)) + (z2 + 2*zeta)*(-(m*(z2 - 2*zeta_conj)) + z2*(z2 - zeta_conj)*zeta_conj))

		coeff2 = (-2*m**2*(z2 - 2*zeta) + 3*m*z2*zeta*(z2 - 2*zeta_conj) + z2**2*zeta*zeta_conj*(-z2 + zeta_conj) + dm*z2*(-2*m - z2*zeta + 2*z2*zeta_conj + 2*zeta*zeta_conj))

		coeff1 = - ((dm - m)*z2*(dm*z2 + m*(z2 - 4*zeta) - z2*zeta*(z2 - 2*zeta_conj)))

		coeff0 = ((dm - m)**2*z2**2*zeta)
		"""

		coeff5 = -(zeta_conj*csum([-z2, zeta_conj]))

		coeff4 = csum([dm*z2, m*csum([z2, -2*zeta_conj]), -csum([2*z2, zeta])*csum([z2, -zeta_conj])*zeta_conj])

		coeff3 = csum([-csum([dm*z2*csum([z2, 2*zeta_conj])]), csum([z2, 2*zeta])*csum([-csum([m*csum([z2, -2*zeta_conj])]), z2*csum([z2, -zeta_conj])*zeta_conj])])

		coeff2 = csum([-2*m**2*csum([z2, -2*zeta]), 3*m*z2*zeta*csum([z2, -2*zeta_conj]), z2**2*zeta*zeta_conj*csum([-z2, zeta_conj]), dm*z2*csum([-2*m, -z2*zeta, 2*z2*zeta_conj, 2*zeta*zeta_conj])])

		coeff1 = -csum([csum([dm, -m])*z2*csum([dm*z2, m*csum([z2, -4*zeta]), -z2*zeta*csum([z2, -2*zeta_conj])])])

		coeff0 = (csum([dm, -m])**2*z2**2*zeta)


	elif (calc == 'caustic'):
		# Specific form of the derived for the planetary caustic frame

		"""
		coeff5 = ((z1 - zeta_conj)*(-z2 + zeta_conj))

		coeff4 = (dm*(-z1 + z2) + m*(z1 + z2 - 2*zeta_conj) + (2*z1 + 2*z2 + zeta)*(z1 - zeta_conj)*(z2 - zeta_conj))

		coeff3 = (-(m*(z1 + z2 + 2*zeta)*(z1 + z2 - 2*zeta_conj)) - (z1**2 + 2*z1*(2*z2 + zeta) + z2*(z2 + 2*zeta))*(z1 - zeta_conj)*(z2 - zeta_conj) + dm*(z1 - z2)*(z1 + z2 + 2*zeta_conj))

		coeff2 = (-2*m**2*(z1 + z2 - 2*zeta) + 3*m*(z1 + z2)*zeta*(z1 + z2 - 2*zeta_conj) + (z2**2*zeta + z1**2*(2*z2 + zeta) + 2*z1*z2*(z2 + 2*zeta))*(z1 - zeta_conj)*(z2 - zeta_conj) + dm*(z1 - z2)*(2*m + z2*zeta + z1*(-2*z2 + zeta - 2*zeta_conj) - 2*z2*zeta_conj - 2*zeta*zeta_conj))

		coeff1 = (-(dm**2*(z1 - z2)**2) + m**2*(z1**2 + 6*z1*z2 + z2**2 - 4*z1*zeta - 4*z2*zeta) + m*(z1*z2*(z2 - 4*zeta) + z1**2*(z2 - zeta) - z2**2*zeta)*(z1 + z2 - 2*zeta_conj) - z1*z2*(2*z2*zeta + z1*(z2 + 2*zeta))*(z1 - zeta_conj)*(z2 - zeta_conj) + dm*(z1 - z2)*(z1**2*(z2 - zeta) - zeta*(4*m + z2*(z2 - 2*zeta_conj)) + z1*(z2**2 - 2*z2*zeta + 2*z2*zeta_conj + 2*zeta*zeta_conj)))

		coeff0 = (dm**2*(z1 - z2)**2*zeta - m**2*(z1 + z2)*(2*z1*z2 - z1*zeta - z2*zeta) - m*z1*z2*(z1*(z2 - zeta) - z2*zeta)*(z1 + z2 - 2*zeta_conj) + z1**2*z2**2*zeta*(z1 - zeta_conj)*(z2 - zeta_conj) - dm*(z1 - z2)*(2*m*(z1*(z2 - zeta) - z2*zeta) + z1*z2*(z1*z2 - z1*zeta - z2*zeta + 2*zeta*zeta_conj)))
		"""

		csum = mm.Utils.complex_fsum

		coeff5 = (csum([z1, -zeta_conj])*csum([-z2, zeta_conj]))

		coeff4 = csum([dm*csum([-z1, z2]), m*csum([z1, z2, -2*zeta_conj]), csum([2*z1, 2*z2, zeta])*csum([z1, -zeta_conj])*csum([z2, -zeta_conj])])

		coeff3 = csum([-(m*csum([z1, z2, 2*zeta])*csum([z1, z2, -2*zeta_conj])), -csum([z1**2, 2*z1*csum([2*z2, zeta]), z2*csum([z2, 2*zeta])])*csum([z1, -zeta_conj])*csum([z2, -zeta_conj]), dm*csum([z1, -z2])*csum([z1, z2, 2*zeta_conj])])

		coeff2 = (-2*m**2*(z1+z2-2*zeta)+3*m*(z1+z2)*zeta*(z1+z2-2*zeta_conj)+(z2**2*zeta+z1**2*(2*z2+zeta)+2*z1*z2*(z2+2*zeta))*(z1-zeta_conj)*(z2-zeta_conj)+dm*(z1-z2)*(2*m+z2*zeta+z1*(-2*z2+zeta-2*zeta_conj)-2*z2*zeta_conj-2*zeta*zeta_conj))

		coeff1 = csum([-(dm**2*csum([z1, -z2])**2), m**2*csum([z1**2, 6*z1*z2, z2**2, -4*z1*zeta, -4*z2*zeta]), m*csum([z1*z2*csum([z2, -4*zeta]), z1**2*csum([z2, -zeta]), -z2**2*zeta])*csum([z1, z2, -2*zeta_conj]), -z1*z2*csum([2*z2*zeta, z1*csum([z2, 2*zeta])])*csum([z1, -zeta_conj])*csum([z2, -zeta_conj]), dm*csum([z1, -z2])*csum([z1**2*csum([z2-zeta]), -zeta*csum([4*m, z2*csum([z2, -2*zeta_conj])]), z1*csum([z2**2, -2*z2*zeta, 2*z2*zeta_conj, 2*zeta*zeta_conj])])])

		coeff0 = csum([dm**2*csum([z1, -z2])**2*zeta, -m**2*csum([z1, z2])*csum([2*z1*z2, -z1*zeta, -z2*zeta]), -m*z1*z2*csum([z1*csum([z2, -zeta]), -z2*zeta])*csum([z1, z2, -2*zeta_conj]), z1**2*z2**2*zeta*csum([z1, -zeta_conj])*csum([z2, -zeta_conj]), -dm*csum([z1, -z2])*csum([2*m*csum([z1*csum([z2, -zeta]), -z2*zeta]), z1*z2*csum([z1*z2, -z1*zeta, -z2*zeta, 2*zeta*zeta_conj])])])


	elif (calc == 'general'):
		# General form of the coefficients
	
		coeff5 = (-zeta_conj + z1)*(zeta_conj- z2)

		coeff4 = (m*z1 + m*z2 + 2.*(z1**2)*z2 + 2.*z1*(z2**2) + dm*(-z1 + z2) + z1*z2*zeta + (zeta_conj**2)*(2.*z1 + 2.*z2 + zeta) - zeta_conj*(2.*m + (z1 + z2)*(2.*z1 + 2.*z2 + zeta)))

		coeff3 = (dm*(z1**2) - m*(z1**2) - 2.*m*z1*z2 - (z1**3)*z2 - dm*(z2**2) - m*(z2**2) - 4.*(z1**2)*(z2**2) - z1*(z2**3) - 2.*m*z1*zeta - 2.*m*z2*zeta - 2.*(z1**2)*z2*zeta - 2.*z1*(z2**2)*zeta - (zeta_conj**2)*((z1**2) + 2.*z1*(2.*z2 + zeta) + z2*(z2 + 2.*zeta)) + zeta_conj*(2.*dm*(z1 - z2) + 2.*m*(z1 + z2 + 2.*zeta) + (z1 + z2)*((z1**2) + 4.*z1*z2 + (z2**2) + 2.*z1*zeta + 2.*z2*zeta)))

		coeff2 = (-2.*(m**2)*(z1 + z2 - 2.*zeta) - 3.*m*(2.*zeta_conj - z1 - z2)*(z1 + z2)*zeta + dm*(z1 - z2)*(2.*m - 2.*z1*z2 + z1*zeta + z2*zeta - 2.*zeta_conj*(z1 + z2 + zeta)) + (zeta_conj - z1)*(zeta_conj - z2)*((z2**2)*zeta + (z1**2)*(2.*z2 + zeta) + 2.*z1*z2*(z2 + 2.*zeta)))

		coeff1 = ((-dm**2)*((z1 - z2)**2) + (m**2)*((z1**2) + 6.*z1*z2 + (z2**2) - 4.*z1*zeta - 4.*z2*zeta) - m*(2.*zeta_conj - z1 - z2)*(z1*z2*(z2 - 4.*zeta) + (z1**2)*(z2 - zeta) - (z2**2)*zeta) - (zeta_conj - z1)*z1*(zeta_conj - z2)*z2*(2.*z2*zeta + z1*(z2 + 2.*zeta)) + dm*(z1 - z2)*(z1*z2*(z2 - 2.*zeta) + (z1**2)*(z2 - zeta) - (4.*m + (z2**2))*zeta + 2.*zeta_conj*(z2*zeta + z1*(z2 + zeta))))

		coeff0 = (-2.*(m**2)*(z1**2)*z2 - 2.*(m**2)*z1*(z2**2) - m*(z1**3)*(z2**2) - m*(z1**2)*(z2**3) + (m**2)*(z1**2)*zeta + (dm**2)*((z1 - z2)**2)*zeta + 2.*(m**2)*z1*z2*zeta + m*(z1**3)*z2*zeta + (m**2)*(z2**2)*zeta + (zeta_conj**2)*(z1**2)*(z2**2)*zeta + 2.*m*(z1**2)*(z2**2)*zeta + m*z1*(z2**3)*zeta + (z1**3)*(z2**3)*zeta - dm*(z1 - z2)*(2.*m + z1*z2)*(z1*(z2 - zeta) - z2*zeta) - zeta_conj*z1*z2*((2.*dm*(z1 - z2) + z1*z2*(z1 + z2))*zeta + m*(-2.*z1*z2 + 2.*z1*zeta + 2.*z2*zeta)))

	else:
		raise ValueError('Not able to retreive coefficients for {} case'.
							 format(calc))

	coefficients = np.array([coeff5, coeff4, coeff3, coeff2, coeff1, coeff0])

	return coefficients






