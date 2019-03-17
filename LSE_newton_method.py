import sys
import matplotlib.pyplot as plt
import numpy as np

def transpose(m,dim):
	trans = []
	for i in range(dim):
		temp = []
		for j in range(len(m)):
			temp.append(m[j][i])
		trans.append(temp)
	return trans

def matrix_mul(m1,m2):
	if len(m1[0]) != len(m2):
		print("matrix_mul dim error!")
		sys.exit(1)

	result = []
	for i in range(len(m1)):
		temp = []
		for j in range(len(m2[0])):
	 		a = 0
		 	for k in range(len(m2)):
		 		a = a + m1[i][k] * m2[k][j]
		 	temp.append(a)
		result.append(temp)

	return result

def get_LU(m):
	U = m
	list_mul = []
	for k in range(len(U)-1):
		for i in range(k+1,len(U)):
			mul = U[i][k]/U[k][k]
			list_mul.append(mul)
			for j in range(len(U[0])):
				U[i][j] = U[i][j] - mul*U[k][j]

	
	index = 0
	L = []
	for i in range(len(m)):
		if i == 0:
			L.append([1])
		else:
			L.append([list_mul[index]])
			index = index + 1

	for i in range(1,len(m)):
		temp = []
		for j in range(len(m)):
			if i > j:
				L[j].append(0)
			elif i == j:
				L[j].append(1)
			else:
				L[j].append(list_mul[index])
				index = index + 1
	return L,U

def get_inv_column(L,U,n):
	b = []

	for i in range(len(L)):
		if i == n:
			b.append(1)
		else:
			b.append(0)

	y = [0 for i in range(len(b))]
	
	for i in range(len(b)):
		mul = 0
		for j in range(len(b)):
			mul = mul + L[i][j] * y[j]

		y[i] = b[i]-mul

	#print(y)
	x = [0 for i in range(len(b))]

	for i in range(len(b)-1,-1,-1):
		mul = 0
		for j in range(len(b)):
			mul = mul + U[i][j] * x[j]

		x[i] = (y[i]-mul)/U[i][i]

	return x

def LSE(data,base,lse_lam):
	A = []
	for i in range(len(data)):
		temp =[]
		for j in range(base-1,-1,-1):
			temp.append(data[i][0] ** j)
		A.append(temp)
	b = []

	for i in range(len(data)):
		b.append([data[i][1]])

	At = transpose(A,base)

	AtA = matrix_mul(At,A)

	for i in range(len(AtA)):
		AtA[i][i] = AtA[i][i] + lse_lam


	L,U = get_LU(AtA)
	invt = []
	for i in range(len(L)):
		invt.append(get_inv_column(L,U,i))

	inv = transpose(invt,len(L))

	ans = matrix_mul(inv,At)

	ans = matrix_mul(ans,b)
	return ans
	

def Newton(data,base):
	A = []

	for i in range(len(data)):
		temp =[]
		for j in range(base-1,-1,-1):
			temp.append(data[i][0] ** j)
		A.append(temp)
	b = []

	for i in range(len(data)):
		b.append([data[i][1]])

	At = transpose(A,base)

	AtA = matrix_mul(At,A)

	x0 = [[1] for i in range(len(AtA))]
	Hessian = [[0 for i in range(len(AtA))] for j in range(len(AtA))]


	for i in range(len(AtA)):
		for j in range(len(AtA[0])):
			Hessian[i][j] = AtA[i][j]

	
	L,U = get_LU(Hessian)
	invt = []
	for i in range(len(L)):
		invt.append(get_inv_column(L,U,i))

	Hessian_inv = transpose(invt,len(L))

	gradient_cal1 = matrix_mul(AtA,x0)
	gradient_cal2 = matrix_mul(At,b)

	gradient = []
	for i in range(len(gradient_cal1)):
		gradient.append([gradient_cal1[i][0] - gradient_cal2[i][0]])

	ressult = matrix_mul(Hessian_inv,gradient)

	ans = []

	for i in range(len(x0)):
		ans.append([x0[i][0] - ressult[i][0]])


	return ans

def cal_error(lse_x,newton_x,data):
	lse_error = 0
	for i in range(len(data)):
		x = data[i][0]
		degree_1 = len(lse_x) - 1
		y = 0
		for j in range(len(lse_x)):
			y = y + lse_x[j][0] * (x ** degree_1)
			degree_1 = degree_1 - 1

		lse_error = lse_error + (y-data[i][1])**2

	newton_error = 0
	for i in range(len(data)):
		x = data[i][0]
		degree_2 = len(newton_x) - 1
		y = 0
		for j in range(len(newton_x)):
			y = y + newton_x[j][0] * (x ** degree_2)
			degree_2 = degree_2 - 1

		newton_error = newton_error + (y-data[i][1])**2

	return lse_error,newton_error

def print_result(lse_x,newton_x,data):
	degree_1 = len(lse_x) - 1
	degree_2 = len(newton_x) - 1

	lse_error,newton_error = cal_error(lse_x,newton_x,data)

	print("LSE:")
	print("Fitting line: ",end = '')
	for i in range(len(lse_x)):
		if degree_1 != 0:
			print("%.15fX^%d + " %(lse_x[i][0],degree_1),end = '')
		else:
			print("%.15f" %(lse_x[i][0]))
		
		degree_1 = degree_1 - 1

	print("Total error: %f" %(lse_error))

	print("Newton:")
	print("Fitting line: ",end = '')
	for i in range(len(newton_x)):
		if degree_2 != 0:
			print("%.15fX^%d + " %(newton_x[i][0],degree_2),end = '')
		else:
			print("%.15f" %(newton_x[i][0]))
		degree_2 = degree_2 - 1

	print("Total error: %f" %(newton_error))

def main():
	argument = sys.argv[1:]
	if len(argument) != 3:
		print("You should input all the arguments,including path to testfile, number of bases n and lambda")
		sys.exit(1)


	path = argument[0]
	base = int(argument[1])
	lse_lam = int(argument[2])
	list = []
	with open(path,'r') as f:
		for x in f:
			list.append(x)
	x = []
	y = []
	min_x = 0
	max_x = 0
	data = []
	for i in range(len(list)):
		data.append([float(list[i].split(',')[0]),float(list[i].split(',')[1])])
		x.append(float(list[i].split(',')[0]))
		y.append(float(list[i].split(',')[1]))

		if i == 0:
			min_x = float(list[i].split(',')[0])
			max_x = float(list[i].split(',')[0])
		else:
			if float(list[i].split(',')[0]) < min_x:
				min_x = float(list[i].split(',')[0])
			if float(list[i].split(',')[0]) > max_x:
				max_x = float(list[i].split(',')[0])

	lse_x = LSE(data,base,lse_lam)

	newton_x = Newton(data,base)

	print_result(lse_x,newton_x,data)


	data_plot_x = np.linspace(min_x-1,max_x+1,int(max_x-min_x))
	degree_1 = len(lse_x) - 1
	degree_2 = len(newton_x) - 1
	y1 = data_plot_x * 0
	y2 = data_plot_x * 0
	for i in range(len(lse_x)):
		y1 = y1 + lse_x[i][0] * (data_plot_x ** degree_1)
		degree_1 = degree_1 - 1

	for i in range(len(newton_x)):
		y2 = y2 + newton_x[i][0] * (data_plot_x ** degree_2)
		degree_2 = degree_2 - 1

	plt.figure()
	plt.subplot(2,1,1)
	plt.scatter(x,y,color='r')
	plt.plot(data_plot_x,y1)

	plt.subplot(2,1,2)
	plt.scatter(x,y,color='r')
	plt.plot(data_plot_x,y2)

	plt.show()

main()