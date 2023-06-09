import numpy as np
import matplotlib.pyplot as plt
import scipy

# set the size of the graph
plt.rcParams["figure.figsize"] = (8, 7)

# import the results from msidh_results.csv

data = np.genfromtxt('msidh_results.csv', delimiter=',', names=True)
print(data)

x = data['settings']
y = data['average_time'] / 60

for y_i in y:
    # print two significant digits
    print(y_i)

coef, _= scipy.optimize.curve_fit(lambda t,a,b: a*np.exp(b*t),  x,  y,  p0=(4, 0.1))
def predict(b):
    return coef[0] * np.exp(coef[1] * b)

pred_x = np.array(list(range(4, 128)))
pred_y = predict(pred_x)

# plot the results
plt.figure()
plt.plot(pred_x, pred_y, label='Fitted curve')
plt.scatter(x, y, label='MSIDH', color='orange', marker='x')

plt.xticks([4,8,16,32,64,128])
plt.xlabel('security parameter (bits)')
plt.ylabel('average time for key exchange (minutes)')
plt.title('MSIDH key exchange time')
plt.legend()

plt.show()

# make a second plot
plt.figure()
pred_x = np.array(list(range(4, 256)))
pred_y = predict(pred_x)

plt.plot(pred_x, pred_y, label='MSIDH prediction')
plt.xticks([4,8,16,32,64,128,192,256])
plt.xlabel('security parameter (bits)')
plt.ylabel('predicted time for key exchange (minutes)')
plt.title('MSIDH key exchange time prediction')
plt.show()

'''
Generation timing data for MSIDH
'''

x = np.array([4, 8, 16, 32, 64, 128])
y = np.array([0.09961, 0.51857, 1.9015, 25.552, 303.69, 2026.062348])

coef, _= scipy.optimize.curve_fit(lambda t,a,b: a*np.exp(b*t),  x,  y,  p0=(4, 0.1))
def predict_gen(b):
    return coef[0] * np.exp(coef[1] * b)

pred_x = np.array(list(range(4, 128)))
pred_y = predict_gen(pred_x)

# plot the results
plt.figure()
plt.plot(pred_x, pred_y / 60, label='Fitted curve')
plt.scatter(x, y / 60, label='Generation time', color='orange', marker='x')

plt.xticks([4,8,16,32,64,128])
plt.xlabel('security parameter (bits)')
plt.ylabel('time for generation (minutes)')
plt.title('MSIDH generation time')
plt.legend()

plt.show()

# make a second plot
plt.figure()
pred_x = np.array(list(range(4, 256)))
pred_y = predict_gen(pred_x)

plt.plot(pred_x, pred_y / 60, label='Prediction')
plt.xticks([4,8,16,32,64,128,192,256])
plt.xlabel('security parameter (bits)')
plt.ylabel('predicted time for generation (minutes)')
plt.title('MSIDH generation time prediction')
plt.show()






