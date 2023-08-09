%%%-------------------------------------------------------------------
%% @doc blas_validation public API
%% @end
%%%-------------------------------------------------------------------

-module(blas_validation_app).

-behaviour(application).

-export([start/2, stop/1]).

start(_StartType, _StartArgs) ->
    Result = blas_validation_sup:start_link(),
    execute_tests:execute(),
    Result.

stop(_State) ->
    ok.

%% internal functions
