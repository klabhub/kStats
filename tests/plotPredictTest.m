classdef plotPredictTest < matlab.unittest.TestCase
    methods (TestClassSetup)
        function addSourceToPath(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root, ...
                IncludingSubfolders=true));
        end
    end

    methods (Test)
        function testCategoricalCombinations(testCase)
            m = testCase.makeCategoricalModel();
            figureHandle = figure('Visible','off');
            testCase.addTeardown(@() close(figureHandle));

            [prediction,low,high,conditions] = lm.plotPredict(m);

            testCase.verifySize(prediction,[6 1]);
            testCase.verifySize(low,[6 1]);
            testCase.verifySize(high,[6 1]);
            testCase.verifySize(conditions,[6 3]);
            testCase.verifyEqual(string(conditions.Properties.VariableNames), ...
                ["A" "B" "S"]);
        end

        function testContinuousValueIsRequired(testCase)
            m = testCase.makeContinuousModel();
            testCase.verifyError(@() lm.plotPredict(m), ...
                'lm:plotPredict:ContinuousValueRequired');
        end

        function testContinuousValuesPairWithCategoricalConditions(testCase)
            m = testCase.makeContinuousModel();
            figureHandle = figure('Visible','off');
            testCase.addTeardown(@() close(figureHandle));

            continuousValues = table((1:3)', 'VariableNames', {'X'});
            [prediction,low,high,conditions] = lm.plotPredict(m, ...
                'continuousValues',continuousValues);

            testCase.verifySize(prediction,[3 1]);
            testCase.verifySize(low,[3 1]);
            testCase.verifySize(high,[3 1]);
            testCase.verifySize(conditions,[3 3]);
            testCase.verifyEqual(conditions.X,(1:3)');
            testCase.verifyEqual(string(conditions.A),["a";"b";"c"]);
        end

        function testSignificantOnlyUsesConfidenceIntervals(testCase)
            m = testCase.makeCategoricalModel();
            figureHandle = figure('Visible','off');
            testCase.addTeardown(@() close(figureHandle));

            [prediction,low,high,conditions] = lm.plotPredict(m, ...
                'sigOnly',true);

            testCase.verifySize(conditions, [numel(prediction) 3]);
            testCase.verifyTrue(all(low > 0 | high < 0));
        end

        function testSignificantDifferenceOnlyUsesIntercept(testCase)
            m = testCase.makeCategoricalModel();
            figureHandle = figure('Visible','off');
            testCase.addTeardown(@() close(figureHandle));

            [~,low,high,conditions] = lm.plotPredict(m, ...
                'sigDiffOnly',true);

            interceptIx = string(m.CoefficientNames) == "(Intercept)";
            interceptLow = m.Coefficients.Lower(interceptIx);
            interceptHigh = m.Coefficients.Upper(interceptIx);
            testCase.verifySize(conditions,[numel(low) 3]);
            testCase.verifyTrue(all(high < interceptLow | low > interceptHigh));
        end

        function testTwoCategoricalPredictionGrid(testCase)
            m = testCase.makeCategoricalModel();
            figureHandle = figure('Visible','off');
            testCase.addTeardown(@() close(figureHandle));

            [prediction,low,high,conditions] = lm.plotPredict(m, ...
                'plot2D',true);

            testCase.verifySize(prediction,[6 1]);
            testCase.verifySize(low,[6 1]);
            testCase.verifySize(high,[6 1]);
            testCase.verifySize(conditions,[6 3]);
            testCase.verifyNumElements(findobj(gca,'Type','image'),1);
            testCase.verifyNumElements(get(gca,'XTick'),2);
            testCase.verifyNumElements(get(gca,'YTick'),3);
        end

        function testTwoCategoricalGridMapsSignificanceToAlpha(testCase)
            m = testCase.makeCategoricalModel();
            figureHandle = figure('Visible','off');
            testCase.addTeardown(@() close(figureHandle));

            lm.plotPredict(m,'plot2D',true,'sigOnly',true);
            imageHandle = findobj(gca,'Type','image');
            alphaData = imageHandle.AlphaData;

            testCase.verifySize(alphaData,[3 2]);
            testCase.verifyFalse(any(isnan(alphaData(:))));
            testCase.verifyGreaterThanOrEqual(alphaData,0.2);
            testCase.verifyLessThanOrEqual(alphaData,1);
        end
    end

    methods (Static, Access=private)
        function m = makeCategoricalModel()
            subjects = categorical(repelem((1:6)',3));
            A = categorical(repmat(["a";"b";"c"],6,1));
            B = categorical(repmat(["x";"y"],9,1));
            Y = 10 + double(A) + .5*double(B) + ...
                repelem([0 .2 -.1 .1 -.2 .1]',3) + ...
                repmat([.1;-.1;.05],6,1);
            T = table(Y,A,B,subjects,'VariableNames', ...
                {'Y','A','B','S'});
            m = fitlme(T,'Y~A+B+(1|S)','DummyVarCoding','reference');
        end

        function m = makeContinuousModel()
            subjects = categorical(repelem((1:6)',3));
            A = categorical(repmat(["a";"b";"c"],6,1));
            X = (1:18)';
            Y = 10 + double(A) + X + ...
                repelem([0 .2 -.1 .1 -.2 .1]',3) + ...
                repmat([.1;-.1;.05],6,1);
            T = table(Y,A,X,subjects,'VariableNames', ...
                {'Y','A','X','S'});
            m = fitlme(T,'Y~A+X+(1|S)','DummyVarCoding','reference');
        end
    end
end
