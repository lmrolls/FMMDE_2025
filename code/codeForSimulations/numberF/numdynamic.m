%% Number of dynamic factors
function qhat = numdynamic(xcmc,dynamic,w)


    %           _______________$$ Hallin Liska criterion $$________________
if strcmp(dynamic,'hallinliska')
    qhat    = hallinliskacriterion2007(xcmc,6,3,w,'p3');
    qhat    = qhat(1);
    if qhat==0
        error('HallinLiska: zero dynamic factors found')
    end
    
	%           _______________$$ Onatski test $$__________________________
elseif strcmp(dynamic,'onatskitest')
    qhat    = onatskitest(xcmc','d');

    %           _______________$$ no criterion $$__________________________
elseif isnumeric(dynamic)
    qhat    = dynamic;
end
%                     _____________________________________________________

end