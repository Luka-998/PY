from sklearn.linear_model import LinearRegression
model = LinearRegression()
import numpy as np

x = np.array([1,2,3,4,5,6])
y = np.array([14,15,16,17,19,20])
x = x.reshape(-1,1)
print(x.shape)
print(y.shape)

print(x)
print(y)

model.fit(x,y)

model1 = LinearRegression().fit(x,y)

