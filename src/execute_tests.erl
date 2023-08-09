-module(execute_tests).
-export([execute/0]).



execute()->
    level1_test(),
    tr_test(),
    ge_test(),
    sp_test(),
    dgetrf_test(),
    ok.
-record(blas, {size, vars, expected}).

build_call_1(#blas{size=S, vars=Call}, Vars)->
    lists:map(
        fun(V)->
            Is_key = orddict:is_key(V, Vars),
            if Is_key ->
                Value = orddict:fetch(V, Vars),
                if 
                    is_list(Value) ->
                        blas:new(chain:ltb(S, Value));
                    (V == n) andalso (S==c orelse S==z) -> 
                        floor(Value / 2);
                    true ->
                        Value
                end;
            true ->
                if V == inc -> 1;
                true -> V end
            end
        end,
        Call
    ).

compare_1(#blas{size=S, vars=Vars}, Lval, Rval)->
    if is_list(Rval)->
        Llist = chain:btl(S, if is_tuple(Lval) -> blas:to_bin(Lval); true-> Lval end),
    
        io:format("Result: l ~p ~p ~p~n", [lists:nth(1, Vars), Llist,  Rval]),
        Llist ==  Rval;
    true->
        io:format("Result: t ~p ~p ~p~n", [lists:nth(1, Vars), Lval,  Rval]),
        Lval == Rval
    end.

run_1(AllVars, Fcts)->
    lists:all(
        fun(Blas=#blas{vars=BlasVars, expected=Expected})->
            Call = build_call_1(Blas, AllVars),
            io:format("Call: ~p~n", [Call]),
            Out = blas:run(list_to_tuple(Call)),
            Result = orddict:from_list(lists:zip([out] ++ BlasVars, [Out] ++ Call)),

            Passed = lists:all(
                fun({Arg, Val})->
                    compare_1(Blas, orddict:fetch(Arg, Result), Val)
                end,
                Expected
            ),

            if not Passed ->
                nok = lists:nth(1, BlasVars);
            true ->
                true
            end
        end,
        Fcts
    ).


level1_test()->
    AllVars = orddict:from_list([
        {n, 4},
        {x, [1,-2,3,-4]},
        {x2, [2,2,2,2]},
        {y, [3,2,0,-2]},
        {a, [1,0]},
        {c, 0.0},
        {s, 1.0}
    ]),

    % Parameters are named in orddict vars. inc is always translated to 1. out is the output of the blas.
    Fcts = [
        % axpy
        #blas{size=s, vars=[saxpy, n, a, x, inc, y, inc], expected=[{y, [4, 0, 3, -6]}]},
        #blas{size=d, vars=[daxpy, n, a, x, inc, y, inc], expected=[{y, [4, 0, 3, -6]}]},
        #blas{size=c, vars=[caxpy, n, a, x, inc, y, inc], expected=[{y, [4, 0, 3, -6]}]},
        #blas{size=z, vars=[zaxpy, n, a, x, inc, y, inc], expected=[{y, [4, 0, 3, -6]}]},
        % copy
        #blas{size=s, vars=[scopy, n, x, inc, y, inc], expected=[{y, [1,-2,3,-4]}]},
        #blas{size=d, vars=[dcopy, n, x, inc, y, inc], expected=[{y, [1,-2,3,-4]}]},
        #blas{size=c, vars=[ccopy, n, x, inc, y, inc], expected=[{y, [1,-2,3,-4]}]},
        #blas{size=z, vars=[zcopy, n, x, inc, y, inc], expected=[{y, [1,-2,3,-4]}]},
        % swap
        #blas{size=s, vars=[sswap, n, x, inc, y, inc], expected=[{x, [3,2,0,-2]}, {y, [1,-2,3,-4]}]},
        #blas{size=d, vars=[dswap, n, x, inc, y, inc], expected=[{x, [3,2,0,-2]}, {y, [1,-2,3,-4]}]},
        #blas{size=c, vars=[cswap, n, x, inc, y, inc], expected=[{x, [3,2,0,-2]}, {y, [1,-2,3,-4]}]},
        #blas{size=z, vars=[zswap, n, x, inc, y, inc], expected=[{x, [3,2,0,-2]}, {y, [1,-2,3,-4]}]},
        % scal
        #blas{size=s, vars=[sscal,  n, a, x, inc], expected=[{x, [1,-2,3,-4]}]},
        #blas{size=d, vars=[dscal,  n, a, x, inc], expected=[{x, [1,-2,3,-4]}]},
        #blas{size=c, vars=[cscal,  n, a, x, inc], expected=[{x, [1,-2,3,-4]}]},
        #blas{size=z, vars=[zscal,  n, a, x, inc], expected=[{x, [1,-2,3,-4]}]},
        #blas{size=c, vars=[csscal, n, a, x, inc], expected=[{x, [1,-2,3,-4]}]},
        #blas{size=z, vars=[zdscal, n, a, x, inc], expected=[{x, [1,-2,3,-4]}]},
        % dot
        %#blas{size=s, vars=[sdot,  n, x, inc, y, inc], expected=[{out, 7}]},
        %#blas{size=d, vars=[ddot,  n, x, inc, y, inc], expected=[{out, 7}]},
        %#blas{size=s, vars=[dsdot, n, x, inc, y, inc], expected=[{out, 7}]},
        %#blas{size=c, vars=[cdotu, n, x, inc, y, inc], expected=[{out, [-1, -10]}]},
        %#blas{size=z, vars=[zdotu, n, x, inc, y, inc], expected=[{out, [-1, -10]}]},
        %#blas{size=c, vars=[cdotc, n, x, inc, y, inc], expected=[{out, [7, 2]}]},
        %#blas{size=z, vars=[zdotc, n, x, inc, y, inc], expected=[{out, [7, 2]}]},
        %#blas{size=s, vars=[sdsdot,n, a, x, inc, y, inc], expected=[{out, 8}]},
        % asum
        %#blas{size=s, vars=[sasum,  n, x, inc], expected=[{out, 10}]},
        %#blas{size=d, vars=[dasum,  n, x, inc], expected=[{out, 10}]},
        %#blas{size=c, vars=[scasum, n, x, inc], expected=[{out, 10}]},
        %#blas{size=z, vars=[dzasum, n, x, inc], expected=[{out, 10}]},
        % sum
        %#blas{size=s, vars=[ssum,  n, x, inc], expected=[{out, -2}]},
        %#blas{size=d, vars=[dsum,  n, x, inc], expected=[{out, -2}]},
        %#blas{size=c, vars=[scsum, n, x, inc], expected=[{out, -2}]},
        %#blas{size=z, vars=[dzsum, n, x, inc], expected=[{out, -2}]}
        % i_min removed
        % i_max removed
        % i_amin
        %#blas{size=s, vars=[isamin,  n, x, inc], expected=[{out, 0}]},
        %#blas{size=d, vars=[idamin,  n, x, inc], expected=[{out, 0}]},
        %#blas{size=c, vars=[icamin, n, x, inc],  expected=[{out, 0}]},
        %#blas{size=z, vars=[izamin, n, x, inc],  expected=[{out, 0}]},
        %i_amax
        %#blas{size=s, vars=[isamax,  n, x, inc], expected=[{out, 3}]},
        %#blas{size=d, vars=[idamax,  n, x, inc], expected=[{out, 3}]},
        %#blas{size=c, vars=[icamax, n, x, inc],  expected=[{out, 1}]},
        %#blas{size=z, vars=[izamax, n, x, inc],  expected=[{out, 1}]},
        %nrm
        %#blas{size=s, vars=[snrm2,  n, x2, inc], expected=[{out, 4}]},
        %#blas{size=d, vars=[dnrm2,  n, x2, inc], expected=[{out, 4}]},
        %#blas{size=c, vars=[scnrm2, n, x2, inc], expected=[{out, 4}]},
        %#blas{size=z, vars=[dznrm2, n, x2, inc], expected=[{out, 4}]},
        %rot
        #blas{size=s, vars=[srot, n, x, inc, y, inc, c, s], expected=[{x, [3,2,0,-2]}, {y, [-1,2,-3,4]}]},
        #blas{size=d, vars=[drot, n, x, inc, y, inc, c, s], expected=[{x, [3,2,0,-2]}, {y, [-1,2,-3,4]}]}
        %#blas{size=c, vars=[csrot, n, x, inc, y, inc, c, s], expected=[{x, [3,2,0,-2]}, {y, [-1,2,-3,4]}]},
        %#blas{size=z, vars=[zdrot, n, x, inc, y, inc, c, s], expected=[{x, [3,2,0,-2]}, {y, [-1,2,-3,4]}]}

    ],

    run_1(AllVars, Fcts).


% file blas3_SUITE

-record(instance, {size, name, expected}).
-define(is_complex(S), (Size==c orelse Size == z)).
-record(blas3, {signature, instances}).

build_call_2(Signature, #instance{size = Size, name = Name}, Vars)->
    Variable_fetcher = fun(Arg_name) ->
        case orddict:find(Arg_name, Vars) of 

            {ok, Value} when is_list(Value) ->
                Split_ratio = if 
                    ?is_complex(Size) andalso (Arg_name == a orelse Arg_name == b orelse Arg_name == c) ->
                        1/2;
                    %?is_complex(Size) ->
                    %    1/2;
                    true ->
                        1
                end,
                Value_list = element(1, lists:split(round(length(Value)*Split_ratio), Value)),
                %io:format("Result: ~p~n", [Value_list]),
                blas:new(chain:ltb(Size, Value_list));

            {ok, Value} when (Arg_name == n) andalso ?is_complex(Size) ->
                floor(Value / 2);

            {ok, Value} -> 
                Value;

            _ when Arg_name == inc->
                1;

            _ ->
                Arg_name
        end
    end,
    [Name] ++ lists:map(Variable_fetcher, Signature).

compare_2({instance, Size, Name, _}, Lval, Rval)->
    Lval_comparable = if
        is_list(Rval) -> chain:btl(Size, blas:to_bin(Lval));
        true          -> Lval
    end,

    io:format("Result: ~p ~p ~p~n", [Name, Lval_comparable, Rval]),
    Lval_comparable == Rval.

run_2(Var_table, #blas3{ signature = Signature, instances = Instances})->
    Test_instance = fun (Instance=#instance{expected = Expected}) ->
        Call = build_call_2(Signature, Instance, Var_table),
        %io:format("Calling ~p~n", [Call]),
        Out  = blas:run(list_to_tuple(Call)),
        State= orddict:from_list(lists:zip([out, name] ++ Signature, [Out] ++ Call)),
        true = lists:all(
            fun({Arg_name, Expected_val})->
                compare_2(Instance, orddict:fetch(Arg_name, State), Expected_val)
            end,
            Expected
        )
    end,

    lists:all(Test_instance, Instances).


ge_test()->
    Vars = orddict:from_list([
        {alpha, [1,0]},
        {beta,  [2,1]},
        {order, blasRowMajor},
        {transA, blasNoTrans},
        {transB, blasTrans},
        {n, 4},
        {a, [1, 3, 2, 4,  5, 3, 9, 2,  4, 1, 2,-2,  2, 2,-1, 1]},
        {b, [2, 2,-1, 1,  2, 2, 1, 5,  7, 2, 5, 0,  1, 3, 2, 4]},
        {c, [1, 0, 0, 0,  0, 1, 0, 0,  0, 0, 1, 0,  0, 0, 0, 1]},

        {x, [1, -4, 2, 2]},
        {y, [1, -4, 1, 2]}
    ]),

    Fcts = [
        % gemm
        #blas3{
            signature = [ order, transA, transB, n, n, n, alpha, a, n, b, n, beta, c, n],
            instances = [
              #instance{size = s, name = sgemm,   expected = [{c, [12, 30, 23, 30, 9, 37, 86, 40, 6, 2, 42, 3, 10, 12, 13, 12]}]},
              #instance{size = d, name = dgemm,   expected = [{c, [12, 30, 23, 30, 9, 37, 86, 40, 6, 2, 42, 3, 10, 12, 13, 12]}]},
              #instance{size = c, name = cgemm,   expected = [{c, [-8, 7, -22, 22, -8, 25, 3, 63]}]},
              #instance{size = z, name = zgemm,   expected = [{c, [-8, 7, -22, 22, -8, 25, 3, 63]}]}
            ]
        },
        %gemv
        #blas3{
            signature = [ order, transA, n, n, alpha, a, n, x, inc, beta, y, inc],
            instances = [
              #instance{size = s, name = sgemv,   expected = [{y, [3.0,7.0,2.0,-2.0]}]},
              #instance{size = d, name = dgemv,   expected = [{y, [3.0,7.0,2.0,-2.0]}]},
              #instance{size = c, name = cgemv,   expected = [{y, [15.0,4.0,31.0,10.0]}]},
              #instance{size = z, name = zgemv,   expected = [{y, [15.0,4.0,31.0,10.0]}]}
            ]
        },
        %ger
        #blas3{
            signature = [ order, n, n, alpha, x, inc, y, inc, a, n],
            instances = [
              #instance{size = s, name = sger,  expected = [{a, [2, -1, 3, 6, 1, 19, 5, -6, 6, -7, 4, 2, 4, -6, 1, 5]}]},
              #instance{size = d, name = dger,  expected = [{a, [2, -1, 3, 6, 1, 19, 5, -6, 6, -7, 4, 2, 4, -6, 1, 5]}]}
              #instance{size = c, name = cgeru, expected = [{a, [-14, -5, 11, 2, 15, -3, 7, 8]}]},
              #instance{size = z, name = zgeru, expected = [{a, [-14, -5, 11, 2, 15, -3, 7, 8]}]},
              #instance{size = c, name = cgerc, expected = [{a, [18, 3, -5, -2, -1, 13, 15, 0]}]},
              #instance{size = z, name = zgerc, expected = [{a, [18, 3, -5, -2, -1, 13, 15, 0]}]}
            ]
        }
    ],

    lists:all(fun (Blas) -> run_2(Vars, Blas) end, Fcts).


sp_test()->
     Vars = orddict:from_list([
        {alpha, [2,0]},
        {beta,  [3,1]},
        {uplo, blasUpper}, 
        {order, blasRowMajor},
        {n, 3},
        {a, [1, 3, 2, 
                3, 9,
                    2]},

        {x, [1, 2, -1]},
        {y, [0, -4, 1]}
     ]),

      Fcts = [
        % spmv
        #blas3{
            signature = [ order, uplo, n, alpha, a, x, inc, beta, y, inc],
            instances = [
              #instance{size = s, name = sspmv,   expected = [{y, [10, -12, 39]}]},
              #instance{size = d, name = dspmv,   expected = [{y, [10, -12, 39]}]}
            ]
        },

        #blas3{
            signature = [ order, uplo, n, alpha, x, inc, a],
            instances = [
              #instance{size = s, name = sspr,   expected = [{a, [3, 7, 0, 11, 5, 4]}]},
              #instance{size = d, name = dspr,   expected = [{a, [3, 7, 0, 11, 5, 4]}]}
            ]
        }, 

         #blas3{
            signature = [ order, uplo, n, alpha, x, inc, y, inc, a],
            instances = [
              #instance{size = s, name = sspr2,   expected = [{a, [1, -5, 4, -29, 21, -2]}]},
              #instance{size = d, name = dspr2,   expected = [{a, [1, -5, 4, -29, 21, -2]}]}
            ]
        }
      ],

      
    lists:all(fun (Blas) -> run_2(Vars, Blas) end, Fcts). 

tr_test()->
     Vars = orddict:from_list([
        {order, blasRowMajor},
        {side, blasLeft},
        {uplo, blasUpper},
        {trans, blasNoTrans},
        {diag, blasNonUnit},
        {n, 4},
        {alpha, [1,-1]},
        {a, [1,0,3,0, 0,5,0,7, 0,0,9,0, 0,0,0,0]},
        {aa, [1,0,0,1, 0,0,1,0]},
        {b, [1,2,3,4, 0,1,6,7, 0,0,1,9, 0,0,0,1]},
        {x, [-1,2,-3,4]},
        {xx, [1,0,1,0]},
        {y, [1,2,3,4]},
        {a_trsm, [1,0,0,0, 0,1,0,0, 0,0,0,0, 0,0,0,1]}
     ]),

      Fcts = [
        #blas3{
            signature = [ order, uplo, trans, diag, n, a, n, x, inc],
            instances = [
                
                % trmv
                #instance{size = s, name = strmv, expected = [{x, [-10,38,-27,0]}]},
                #instance{size = d, name = dtrmv, expected = [{x, [-10,38,-27,0]}]},
                #instance{size = c, name = ctrmv, expected = [{x, [-10,14,-28,-21]}]},
                #instance{size = z, name = ztrmv, expected = [{x, [-10,14,-28,-21]}]}
            ]
        },
        #blas3{
            signature = [ order, uplo, trans, diag, n, b, n, y, inc],
            instances = [
                 % trsv
                #instance{size = s, name = strsv, expected = [{y, [-260, 172, -33, 4]}]},
                #instance{size = d, name = dtrsv, expected = [{y, [-260, 172, -33, 4]}]}
            ]
        },
         #blas3{
            signature = [ order, uplo, trans, diag, n, aa, n, xx, inc],
            instances = [
                 % trsv
                #instance{size = c, name = ctrsv, expected = [{xx, [1,-1,1,0]}]},
                #instance{size = z, name = ztrsv, expected = [{xx, [1,-1,1,0]}]}
            ]
        }
      ],
      
    lists:all(fun (Blas) -> run_2(Vars, Blas) end, Fcts). 


dgetrf_test()->
    A   = blas:new(float64, [1, 2, -3, 1, 2, 4, 0, 7, -1, 3, 2, 0]),
    IPV = blas:new(int64, [0,0,0]),

    ok  = blas:run({dgetrf, blasRowMajor, 3, 4, A, 4, IPV}),

    [2,3,3] = blas:btl(int32, blas:to_bin(IPV)),
    
    A_expected = [2.0,4.0,0.0,7.0,-0.5,5.0,2.0,3.5,0.5,0.0,-3.0,-2.5],
    A_expected =:= blas:btl(float64, blas:to_bin(A)).

