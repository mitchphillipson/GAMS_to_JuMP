# GAMS_to_JuMP
 A collection of models written in both GAMS and Julia JuMP


01 - Introduction
- Linear climate model


07 - PE Models
- marshall - MCP model formulated as both MCP and NLP
- qpdispatch - MCP model formulated as both MCP and NLP

12 - Expected Utility
- insurance - MCP model formulated as both MCP and NLP
- NashInsurance - NLP model. Note: Found bug in GAMS code while debugging this. In GAMS there is a macro defined as
    `$macro U(C)	(C**rho/rho)`
    however, this gives an unexpected result when run as
    `U(1+L)`
    the expected result is 
    $$\frac{(1+L)^\rho}{\rho}$$ 
    however GAMS expands the macro as
    $$1 + \frac{L^\rho}{\rho}$$
    which is an issue. The solution is to change the macro to 
    `$macro U(C)	((C)**rho/rho)` to force the parentheses.

    The solutions don't match perfectly. There are about 40 solve iterations and the middle solutions don't match.

13 - Mor Expected Utility
- NashInsurance - NLP model. Similar to the model in 12, still debugging. Iterations aren't matching GAMS.

22 - NonConvexity
- game22 - MIP model of triangle sum game
- multiple - NLP plotting a risky route