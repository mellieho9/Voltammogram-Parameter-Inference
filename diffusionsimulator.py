## Populate number of cells based on user input
print("Enter the number of cells you want to start with")
num_cells = input()
cont_at_all_cells=[]
for x in num_cells:
	cont_at_all_cells[x,0]=0.5

## Obtains the parameters to begin simulating diffusion
print("Enter the rate of diffusion and diffusion coefficient")
R = input()
D = input()

## Simulates diffusion
t = 1
while True:
	cont_at_all_cells[0,t] = cont_at_all_cells[x,t-1] + R*t + D*(cont_at_all_cells[1,t-1] - cont_at_all_cells[0,t-1])*t
	for y in (len(cont_at_all_cells) - 2):
		cont_at_all_cells[y+1,t] = cont_at_all_cells[y+1,t-1] + D*(cont_at_all_cells[y+2,t-1] - cont_at_all_cells[y+1,t-1])*t + D*(cont_at_all_cells[y+1,t-1] - cont_at_all_cells[y,t-1])*t
	cont_at_all_cells[-1,t] = cont_at_all_cells[-1,t-1] + D*(cont_at_all_cells[-2,t-1] - cont_at_all_cells[-1,t-1])*t
	t+=1
	if cont_at_all_cells[0,t]==cont_at_all_cells[-1,t]:
		break

## Print results
print("The time it takes for all cells to completely diffuse is: %f seconds",t)
print("The resulting concentration is: %f", cont_at_all_cells[-1,t])
