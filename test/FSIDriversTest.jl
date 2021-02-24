module FSIDriversTest

using GridapFSI.Drivers
using Test

problem = FSIProblem{:not_implemented_problem}()
@test_throws ErrorException execute(problem)
kwargs = Dict(:a=>1.0)
@test Drivers._get_kwarg(:a,kwargs) == 1.0
@test Drivers._get_kwarg(:b,kwargs,2.0) == 2.0
@test_throws UndefVarError Drivers._get_kwarg(:b,kwargs)

end
