# E2CoPre
This is the code repository for the paper $E^2CoPre$: Energy Efficient and Cooperative Collision Avoidance for UAV Swarms with Trajectory Prediction. 

Acknowledgement: this code is partially build on previous works as follows. 
- Chenyang Xu and Jerry L. Prince, "Gradient Vector Flow: A New external force for Snakes, for Gradient Vector Flow. 
- D.Kroon University of Twente (July 2010), for Snake energy optimizaiton.  

Parameter Setting: 
Our simulation set-up is based on the application scenario illustrated in Fig.~\ref{fig:scenario}. Obstacles intrude between any two waypoints on the static path of a swarm. The two waypoints are the start and target positions of the swarm during avoidance.  As illustrated in the previous sections, the trajectories of the UAVs only depend on the relative positions of the UAVs and the obstacles. Therefore, the simulation environment is orientation invariant. For simplicity, the simulation setup is designed as follows. 
The simulation environment is a $300\times300$ square. The UAV swarm is always spawned at the left side of the square and flies toward the right side. 
The spawn positions of the UAV swarm are the same for each run. 
Without loss of generality, we test two specific scenarios: \textit{Obstacle in Front} where the obstacles fly directly toward the UAV swarm and \textit{Obstacle on Side} where the obstacles approach the UAV swarm from the left or right side. 
The two scenarios are illustrated in Fig. \ref{exp-scenario}. 
\begin{figure}[t]
	\centering
	\subfloat[]{\includegraphics[width=0.5\linewidth]{figures/exp-scenario1.pdf}%
		\label{exp_scenario1}}
	\hfil
	\subfloat[]{\includegraphics[width=0.5\linewidth]{figures/exp-scenario2.pdf}%
		\label{exp_scenario2}}
	\caption{Experiment scenarios. (a): \textit{Obstacle in Front}. (b): \textit{Obstacle on Side}. }
 \label{exp-scenario}
\end{figure}

Essential parameters are defined as follows. The initial distance between the swarm and the obstacle is 200 $m$. 
Because the maximum ground speed of industry-level UAVs is between 10 $m/s$ and 17 $m/s$ \cite{DJIIndust}, and the maximum ground speed of consumer-level UAVs is between 5 $m/s$ to 14 $m/s$ \cite{DJIConsum}, we set the swarm speed $v_s=10\ m/s$ in all simulations. The weight of the UAV is set to $1\ Kg$. 
In the simulations, we are interested in the relation between $E^2CoPre$'s performances and the obstacle's velocity rather than extreme tests on the obstacle's velocity. For the sake of simulation, we assume the obstacles are adversarial UAVs that have the same velocity range with the swarm. The velocity of obstacle $v_{obs}$ varies from 1 to 10 $m/s$. The Lidar sensor's sensing range is set to $100\ m$. An obstacle is detected by a UAV when the distance between them is smaller than $100\ m$. The obstacle detected by any UAVs will be known to the whole swarm through information exchange. Each run of the simulation starts when the swarm and the obstacles are spawn and finishes when the distances between all UAVs arrive at their target positions. 
During each run, the collision avoidance starts when the distance between the obstacle and any UAV is smaller than the avoiding distance $50\ m$. 
The length of one planning step is set to $|S| = v_s\times1s$. We predict for 10 planning steps in trajectory planning which is $100\ m$, as it's the length of sensing range of UAVs. 
The radius of the \textit{Protection Bubble} of an obstacle is set to $d_{safe} = 20\ m$. 
The threshold distance for V2O collisions is set to $d_{obs}=d_{safe}-|S|=10/ m$, as the \textit{Protection Bubble} guarantees that UAVs get closer to obstacles with no more than one planning step. 
Lastly, the threshold distance for V2V collisions is set to $d_{v2v}=5\ m$. 
For swarm formation, the swarm members are equally spaced in a circle with a radius $\tau=20\ m$. The swarm size $N$ varies from 2 to 10. The resolution of the simulations is 1 $s$. 
