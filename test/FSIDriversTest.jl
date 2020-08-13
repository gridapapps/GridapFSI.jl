module FSIDriversTest

using GridapFSI.FSIDrivers
using Test

problem = Problem{:not_implemented_problem}()
@test_throws ErrorException execute(problem)
kwargs = Dict(:a=>1.0)
@test FSIDrivers._get_kwarg(:a,kwargs) == 1.0
@test FSIDrivers._get_kwarg(:b,kwargs,2.0) == 2.0

end
