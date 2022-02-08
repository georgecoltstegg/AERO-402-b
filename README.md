# Cislunar Navigation and Communication Constellation (CNCC)

AERO 402 Team B Github Repository

![Alt text](/Pictures/MissionDescription.png)

# Steps to run the STK Scenario Generation:

1. Install STK 12.2.0 from the AGI website.
2. Install the provided STK license from Dr. Selva on Canvas.
3. Using your python cmd, pip install the python API package from the STK install folder under ~\AGI\STK 12\bin\AgPythonAPI
4. Run STK and create a new scenario.
5. Run the python code CNCC_ConstellationDesignCode.py from this github.
6. Watch the scenario be populated with our constellation! 

![Alt text](/Pictures/STK_initialDesign.png)

# Simulator Design 
## Orbital Propagator
This role is directly necessary for all of the requirements. It is required to evaluate the dilution of precision that establishes position accuracy through system geometry. It will assist in determining data throughput by evaluating satellite position over time such that the link budget can be calculated, and the data rate can be evaluated. The orbital propagator is necessary to evaluate the number of required satellites and to evaluate the amount of stationkeeping required by each satellite over the mission duration, which plays into cost.

## Link Budget Calculator
This role is directly necessary for high-level requirements. Dilution of precision directly determines position accuracy. The link budget equation determines data throughput. The link budget equation also tells us the required power for simultaneous connections. The required power and the antenna types and configurations play a role in the budget. It will also assist in determining orbit patterns which have an enormous effect on the budget.

## Architecture Evaluator
This role primarily serves to drastically ease the effort associated with making major changes to the design by automatically handling the consequences of the design choices. It will play a key role in allowing the optimizer to make design decisions and iterate toward the best solution, and will also allow the team to quickly test designs they personally hypothesize.

## Customer and Network Cost Modeler
Allows us to ensure we fulfill the high-level requirements in relation to budget. It also acts as a penalty when designing different systems with the optimizer that will allow us to choose the best possible design.

## Event Handler
This role is responsible for determining WHEN and HOW certain events take place. The nature of this role is fundamentally arbitrary, but should give us satisfaction that the high level requirements are being met.

