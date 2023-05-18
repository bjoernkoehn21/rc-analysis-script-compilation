This document lists some syntax and style guidelines that should be followed in the scripts. The guidelines heavily lean on  [Matlab style guideline 2.0](https://uk.mathworks.com/matlabcentral/fileexchange/46056-matlab-style-guidelines-2-0).

## Naming Constants, Variables and Functions
In general the name should reflect what the variable, constant or functions in question is used for. In particular uninformative counter names such as 'i' or 'j' should be avoided as should be names like 'x' or 'y' for elements in a calculation - if not used in general purpose functions or as scratch variables (who are not used outside a few lines of code).
Separation of words by underscores is not recommended. For explanation, see [Matlab style guideline 2.0](https://uk.mathworks.com/matlabcentral/fileexchange/46056-matlab-style-guidelines-2-0). So, where not misleading, the following naming convention should be followed:
*  constant should be named in `UPPER_CASES` with words separated by an underscore, and
*  variables should be named in `camelCase`,
*  functions should be named in `camelCase`,
*  structures should be named in `camelCase`.

Constants should be declaired first in a function, separated by an empty line from the help or following code.

Here is a code example demonstrating the application of the naming conventions:
```matlab
function costOutput = householdCostCalculation(foodCost, rentCost, nDays)
% costOutput = householdCostCalculation(foodCost, rentCost, nDays)
% Calculates the total household costs for a certain timespan.
%
% INPUT:
%   foodCost - daily costs incurred for food
%   rentCost - dayli costs related to rent
%   nDays - number of days of timespan of interest
% OUTPUT:
%   costOutput  - struct containing days and associatied costs
%     costOutput.day -  date of the period in question
%     costOutput.cost - household costs from now to the day in question

SAFETY_BUFFER = 30; % This is a constant

costOutput.day = datetime('now') + days(0:(nDays-1)); % This is a struct, calculated using a function

dailyCost = foodCost + rentCost;
costOutput.cost = dailyCost * (1:nDays)
costOutput.cost(end) += SAFETY_BUFFER;
end
```

## General coding rules
* All lines of code
 * **should** be no longer than `100` characters.
 * **must** be no longer than `120` characters.
* Insert whitespaces
 * after `,`, i.e. `a = [1, 2, 3];`
 * after `;` if part of a vector or matrix, i.e. `a = [1; 2; 3];`, but not at the end of the line.
 * before and after `=`, i.e., `a = 1`
 * before and after `==`, i.e., `find(a == b)`
* Insert no whitespaces before and after
 * `(` and `)`
 * `{` and `}`
 * `[` and `]`
 * `;` if not part of a vector or matrix
 * `:`
 * `.`
* No line breaks
 * between function name and opening bracket, i.e., a line break should look like this: `result = myFunction(...`
* Files **must** always end with an empty line. The motivation for that can be found [here](http://stackoverflow.com/questions/2287967/why-is-it-recommended-to-have-empty-line-in-the-end-of-file).

## Rules on commented lines
 * No empty line after `%%` or `%`, i.e., comments are never followed by an empty line, no matter what follows (except when writing help for a function, see [the following section](#adding-help-to-a-matlab-function)).
 * Always insert an empty line between actual code and following`%` or `%%`.
 * Text in a comment should start with a capital letter.
 * Text in the second line of a multi-line comment starts with a lowercase letter.
``` matlab
            %% Definition of constants
            CONSTANT_ONE = 1;
            CONSTANT_ONE = 1e-0;

            %% Paths to extra toolboxes and my own functions
            % This is a very long comment that is indeed so long that it takes up so much space that
            % it runs over two lines.
            addpath here/follows/a/path/to/a/file;
            addpath(genpath('here/follows/a/path/to/a/folder'));

            % Some small comment explaining the next block
            N_RANDOM_NETWORKS = 2500;
            N_NODES = 50;
```

## Adding help to a MATLAB function
If a comment is added right after the function definition and closed with an empty line, this whole comment will be returned when `help` is called for this specific function.
The first line of the comment should repeat the interface of the function.
An example of possible help taken form [here](http://de.mathworks.com/help/matlab/matlab_prog/add-help-for-your-program.html) is given below.
You can also find further help concerning function documentation there.
Note the empty line between the help and the function's code.
Also note the indentation shown in the example below:
In the function documentation every subsection should be two whitespaces (`  `) farther to the right, than the section it is derived from.

```matlab
function c = addMe(a, b)
%   C = ADDME(A) adds A to itself.
%   C = ADDME(A, B) adds A and B together.
%
% ADDME  Add two values together.
%   INPUT:
%     a - first summand of the
%       addition performed
%     b - second summand
%   See also SUM, PLUS.

switch nargin
    case 2
        c = a + b;
    case 1
        c = a + a;
    otherwise
        c = 0;
end
```
