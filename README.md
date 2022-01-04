# Voltammogram-Parameter-Inference
A neural network model that infers parameters for glucose fuel cells to generate electricity current

## Background
Glucose is one of the vital energy sources for living species. It is formed after carbohydrates are broken down and are generally oxidized to CO2 through various aerobic metabolic pathways. As it is abundant in the human body, scientists have utilized it in fuel cells to power medical implants i.e. pacemakers or brain stimulators. 

A glucose fuel cell generates electricity through the oxidation of glucose and reduction of oxygen. Unlike batteries, which need to be replaced once their stored chemicals are exhausted, fuel cells work as long as the fuel is replenished. Fuel cells that run on glucose could be a good alternative to lithium-ion batteries for powering implanted devices. 

Upon experimentation, electrical signals produced from fuel cells are highly obscure, raising questions on whether it will provide sufficient electricity to run one device. This project seeks to deconvolute those instances by:
1. Simulate a cyclic voltammogram - an electrochemical device that measures the current by cycling the potential of a working electrode. 
2. Find a gradient descent through the use of real data on the simulated cyclic voltammogram.
3. Use variants of NeuralODE to produce an inverse model that compute the necessary properties that a glucose fuel cell ought to have to produce sufficient electricity.

## Progress
To be updated

## License
[MIT](https://choosealicense.com/licenses/MIT/)
