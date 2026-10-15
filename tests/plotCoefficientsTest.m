classdef plotCoefficientsTest < matlab.unittest.TestCase
    % Tests for represented coefficient totals in lm.plotCoefficients.

    methods (TestClassSetup)
        function addSourceToPath(testCase)
            root = fileparts(fileparts(mfilename('fullpath')));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(root, ...
                IncludingSubfolders=true));
        end
    end

    methods (Test)
        function testSingleCategoricalReference(testCase)
            m = testCase.makeModel("reference",1);
            [estimate,low,high,labels] = testCase.plot(m);

            coefficientNames = string(m.CoefficientNames);
            beta = m.Coefficients.Estimate;
            represented = coefficientNames ~= "(Intercept)";
            expected = beta(represented);

            testCase.verifyEqual(labels, coefficientNames(represented).');
            testCase.verifyEqual(estimate, expected, AbsTol=1e-10);
            testCase.verifyTrue(all(low <= estimate));
            testCase.verifyTrue(all(high >= estimate));

            target = coefficientNames == "A_b";
            contrast = double(coefficientNames == "A_b");
            se = sqrt(contrast*m.CoefficientCovariance*contrast');
            critical = (m.Coefficients.Upper(target) - ...
                m.Coefficients.Estimate(target))/m.Coefficients.SE(target);
            targetIndex = labels == "A_b";
            testCase.verifyEqual(low(targetIndex), ...
                estimate(targetIndex) - critical*se, AbsTol=1e-10);
            testCase.verifyEqual(high(targetIndex), ...
                estimate(targetIndex) + critical*se, AbsTol=1e-10);
        end

        function testSingleCategoricalEffects(testCase)
            m = testCase.makeModel("effects",1);
            [estimate,low,high,labels] = testCase.plot(m);

            coefficientNames = string(m.CoefficientNames);
            beta = m.Coefficients.Estimate;
            represented = coefficientNames( ...
                startsWith(coefficientNames,"A_"));
            expected = [beta(ismember(coefficientNames,represented)); ...
                -sum(beta(ismember(coefficientNames,represented)))];

            testCase.verifyEqual(labels,[represented(:); "A_c"]);
            testCase.verifyEqual(estimate, expected, AbsTol=1e-10);
            testCase.verifyEmpty(intersect(labels,"(Intercept)"));
            testCase.verifyTrue(all(low <= estimate));
            testCase.verifyTrue(all(high >= estimate));
        end

        function testOmittedEffectsLevel(testCase)
            m = testCase.makeModel("effects",1);
            [estimate,low,high,labels] = testCase.plot(m,"A_c");

            coefficientNames = string(m.CoefficientNames);
            beta = m.Coefficients.Estimate;
            contrast = -double(startsWith(coefficientNames,"A_"));
            expected = contrast * beta;

            testCase.verifyEqual(labels,"A_c");
            testCase.verifyEqual(estimate,expected,AbsTol=1e-10);
            testCase.verifyLessThanOrEqual(low,estimate);
            testCase.verifyGreaterThanOrEqual(high,estimate);
        end

        function testOmittedEffectsLevelWithoutIntercept(testCase)
            m = testCase.makeModel("effects",1);
            [estimate,~,~,labels] = testCase.plotWithoutIntercept(m,"A_c");

            coefficientNames = string(m.CoefficientNames);
            beta = m.Coefficients.Estimate;
            expected = -sum(beta(startsWith(coefficientNames,"A_")));

            testCase.verifyEqual(labels,"A_c");
            testCase.verifyEqual(estimate,expected,AbsTol=1e-10);
        end

        function testTwoCategoricalsReference(testCase)
            m = testCase.makeModel("reference",2);
            [estimate,low,high,labels] = testCase.plot(m);

            coefficientNames = string(m.CoefficientNames);
            beta = m.Coefficients.Estimate;
            represented = coefficientNames ~= "(Intercept)";
            expected = beta(represented);

            testCase.verifyEqual(labels,coefficientNames(represented).');
            testCase.verifyEqual(estimate, expected, AbsTol=1e-10);
            testCase.verifyTrue(all(low <= estimate));
            testCase.verifyTrue(all(high >= estimate));
        end

        function testInteractionTerm(testCase)
            m = testCase.makeInteractionModel("reference");
            coefficientNames = string(m.CoefficientNames);
            interaction = coefficientNames == "A_b:B_y";
            beta = m.Coefficients.Estimate(interaction);

            [estimate,low,high,labels] = testCase.plot(m,"B_y:A_b");

            testCase.verifyEqual(labels,"B_y:A_b");
            testCase.verifyEqual(estimate,beta,AbsTol=1e-10);
            testCase.verifyLessThanOrEqual(low,estimate);
            testCase.verifyGreaterThanOrEqual(high,estimate);
        end

        function testOmittedEffectsInteractionTerm(testCase)
            m = testCase.makeInteractionModel("effects");
            [estimate,low,high,labels] = testCase.plot(m,"A_b:B_y");

            coefficientNames = string(m.CoefficientNames);
            beta = m.Coefficients.Estimate;
            contrast = double(coefficientNames == "A_a:B_x");
            expected = contrast * beta;

            testCase.verifyEqual(labels,"A_b:B_y");
            testCase.verifyEqual(estimate,expected,AbsTol=1e-10);
            testCase.verifyLessThanOrEqual(low,estimate);
            testCase.verifyGreaterThanOrEqual(high,estimate);
        end

        function testPartlyOmittedEffectsInteractionTerm(testCase)
            m = testCase.makeInteractionModel("effects");
            [estimate,low,high,labels] = testCase.plot(m,"A_a:B_y");

            coefficientNames = string(m.CoefficientNames);
            beta = m.Coefficients.Estimate;
            contrast = -double(coefficientNames == "A_a:B_x");
            expected = contrast * beta;

            testCase.verifyEqual(labels,"A_a:B_y");
            testCase.verifyEqual(estimate,expected,AbsTol=1e-10);
            testCase.verifyLessThanOrEqual(low,estimate);
            testCase.verifyGreaterThanOrEqual(high,estimate);
        end

        function testOmittedEffectsLevelWithContinuousInteraction(testCase)
            m = testCase.makeContinuousInteractionModel;
            [estimate,low,high,labels] = testCase.plot(m,"A_b:X");

            coefficientNames = string(m.CoefficientNames);
            beta = m.Coefficients.Estimate;
            contrast = -double(coefficientNames == "A_a:X");
            expected = contrast * beta;

            testCase.verifyEqual(labels,"A_b:X");
            testCase.verifyEqual(estimate,expected,AbsTol=1e-10);
            testCase.verifyLessThanOrEqual(low,estimate);
            testCase.verifyGreaterThanOrEqual(high,estimate);
        end

        function testWildcardIncludesOmittedInteractionTerms(testCase)
            m = testCase.makeInteractionModel("effects");

            [~,~,~,labels] = testCase.plot(m);

            coefficientNames = string(m.CoefficientNames);
            coefficientNames = coefficientNames(:);
            expected = [coefficientNames(coefficientNames ~= "(Intercept)"); ...
                "A_b"; "B_y"; "A_a:B_y"; "A_b:B_x"; "A_b:B_y"];
            testCase.verifyEqual(labels,expected);
        end

        function testWildcardMinusExcludesIntercept(testCase)
            m = testCase.makeModel("effects",1);
            [~,~,~,labels] = testCase.plotWithoutIntercept(m,"*-");

            testCase.verifyFalse(any(labels == "(Intercept)"));
            testCase.verifyEqual(labels, ["A_a"; "A_b"; "A_c"]);
        end
    end

    methods (Static, Access=private)
        function m = makeModel(coding,numCategoricals)
            subjects = categorical(repelem((1:6)',3));
            A = categorical(repmat(["a"; "b"; "c"],6,1));
            B = categorical(repmat(["x"; "y"],9,1));
            subjectEffect = repelem([0 .2 -.1 .1 -.2 .1]',3);
            noise = repmat([.1; -.1; .05],6,1);
            Y = 10 + (double(A)-1) + .5*(double(B)-1) + ...
                subjectEffect + noise;
            T = table(Y,A,B,subjects,'VariableNames', ...
                {'Y','A','B','S'});
            if numCategoricals == 1
                m = fitlme(T,'Y~A+(1|S)','DummyVarCoding',coding);
            else
                m = fitlme(T,'Y~A+B+(1|S)','DummyVarCoding',coding);
            end
        end

        function m = makeInteractionModel(coding)
            if nargin < 1
                coding = "reference";
            end
            subjects = categorical(repelem((1:6)',4));
            A = categorical(repmat(["a"; "a"; "b"; "b"],6,1));
            B = categorical(repmat(["x"; "y"; "x"; "y"],6,1));
            subjectEffect = repelem([0 .2 -.1 .1 -.2 .1]',4);
            noise = repmat([.1; -.1; .05; -.05],6,1);
            Y = 10 + (double(A)-1) + .5*(double(B)-1) + ...
                2*(double(A)-1).*(double(B)-1) + subjectEffect + noise;
            T = table(Y,A,B,subjects,'VariableNames', ...
                {'Y','A','B','S'});
            m = fitlme(T,'Y~A*B+(1|S)','DummyVarCoding',coding);
        end

        function m = makeContinuousInteractionModel
            subjects = categorical(repelem((1:6)',4));
            A = categorical(repmat(["a"; "a"; "b"; "b"],6,1));
            X = repmat([-1; -.3; .4; 1],6,1);
            subjectEffect = repelem([0 .2 -.1 .1 -.2 .1]',4);
            noise = repmat([.1; -.1; .05; -.05],6,1);
            Y = 10 + 2*X + (double(A)-1) + ...
                1.5*(double(A)-1).*X + subjectEffect + noise;
            T = table(Y,A,X,subjects,'VariableNames',{'Y','A','X','S'});
            m = fitlme(T,'Y~A*X+(1|S)','DummyVarCoding','effects');
        end
    end

    methods (Static, Access=private)
        function [estimate,low,high,labels] = plot(m,fe)
            if nargin < 2
                fe = "*-";
            end
            figureHandle = figure('Visible','off');
            cleanup = onCleanup(@() close(figureHandle));
            [estimate,low,high] = lm.plotCoefficients(m,fe,[],"none");
            labels = string(get(gca,'XTickLabel'));
            labels = labels(:);
        end

        function [estimate,low,high,labels] = plotWithoutIntercept(m,fe)
            figureHandle = figure('Visible','off');
            cleanup = onCleanup(@() close(figureHandle));
            [estimate,low,high] = lm.plotCoefficients(m,fe,[],"none");
            labels = string(get(gca,'XTickLabel'));
            labels = labels(:);
        end
    end
end
