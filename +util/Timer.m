classdef Timer < matlab.mixin.Copyable  % Inherits copy functionality
    %{
The Timer class provides a flexible and convenient way to measure time in your MATLAB programs. It offers functionality similar to the built-in tic and toc functions but with additional features for tracking multiple time intervals ("lapses") and generating detailed reports.

Key Features:

1. Start and Stop: The start() and stop() methods initiate and terminate the timer, respectively.
2. Lapses: The lapse() method allows you to record intermediate time points within a running timer. This is useful for tracking the duration of different sections of code or events within a process.
3. Duration and Intervals: The duration property provides the total elapsed time, while the intervals property gives the durations between consecutive lapses. Both are formatted as tables for easy readability.
4. Lags: The lags property provides the durations of each lapse relative to the start time.
5. Reporting: The report() method generates formatted output summarizing the elapsed time, intervals, or lags.
6. Copyable: The class inherits from matlab.mixin.Copyable, allowing you to create copies of Timer objects.
    %}
    % This class provides a simple timer functionality, similar to tic/toc, but
    % with the ability to record multiple lapses and control output verbosity.
    % MOz, Feb 25
    properties (SetAccess = protected) % Properties that can be set only within the class
        t0 % Start time of the timer
        tF  % Stop time of the timer
        tN (1,:) datetime = datetime.empty() % Array to store lapse times
    end

    properties (Dependent) % Dependent properties (calculated on the fly)
        duration % Formatted duration string
        lags % lapse durations relative to t0
        intervals % lapse durations relative to tN-1
        n_lapse % Number of lapses recorded

    end

    properties (Dependent, Access = protected) % Dependent, accessible only within the class and its subclasses
        raw_duration % Raw duration as a duration object
        raw_intervals % Raw interval durations as an array of duration objects
        raw_lags % Raw lag durations as an array of duration objects
    end

    properties (Access = protected) % Private properties, accessible only within the class and its subclasses
        isInitiated = false % Flag indicating if the timer has started
        isTerminated = false % Flag indicating if the timer has stopped
    end

    properties(Access = protected, Constant) % Constant properties, accessible only within the class and its subclasses
        duration_types = ["days", "hours", "minutes", "seconds", "milliseconds"]
    end

    methods % Methods of the Timer class
        function self = Timer() % Constructor
            % Creates a Timer object.
            %
            % Usage:
            %   timerObj = Timer();
        end

        function varargout = start(self, varargin) % Start the timer
            % Starts the timer. Optionally displays a message.
            %
            % Usage:
            %   timerObj = start(timerObj); % Start the timer
            %   timerObj = start(timerObj, 'Starting the timer...'); % Start and display a message

            if self.isInitiated && ~self.isTerminated
                error("The existing timer is being reset!");
            end

            self.t0 = datetime("now"); % Record the current time
            self.isInitiated = true; % Set initiation flag
            self.isTerminated = false; % Reset termination flag

            if nargout == 0, varargout = {};
            elseif nargout == 1, varargout{1} = self;
            else, error("Too many output arguments!")
            end

            if nargin > 1, fprintf(varargin{:}); end
        end

        function varargout = stop(self, varargin) % Stop the timer

            if self.isInitiated && ~self.isTerminated % Check if timer was started
                self.tF = datetime("now"); % Record the stop time
                self.isTerminated = true; % Set termination flag
            else
                warning("No timers were initiated. Ignoring.") % Warning if timer wasn't started
                return;
            end

            if nargout == 0, varargout = {};
            elseif nargout == 1, varargout{1} = self;
            else, error("Too many output arguments!")
            end

            if nargin > 1, fprintf(varargin{:}); end
        end

        function varargout = reset(self, varargin)

            if self.isInitiated && ~self.isTerminated

                self.stop();
                self.start(varargin{:});

            else
                warning("No timer was set to reset. Ignoring...");
                return
            end

            if nargout == 0, varargout = {};
            elseif nargout == 1, varargout{1} = self;
            else, error("Too many output arguments!")
            end

            if nargin > 1, fprintf(varargin{:}); end

        end

        function varargout = lapse(self, varargin) % Record a lapse time
            if nargin == 1, varargin{1} = ''; end
            if self.isInitiated || ~self.isTerminated % Check if timer is running
                self.tN(end+1) = datetime("now"); % Record lapse time
            else
                warning("No timers were initiated. Ignoring.") % Warning if timer wasn't started
                varargin = {}; % prevent display
            end

            if nargout == 0, varargout = {};
            elseif nargout == 1, varargout{1} = self;
            else, error("Too many output arguments!")
            end

            if nargin > 1, fprintf(varargin{:}); end
        end

        function d = get.raw_duration(self) % Get raw duration
            % Calculates the raw duration as a duration object.

            d = self.tF - self.t0; % Duration = stop time - start time

        end

        function d = get.duration(self) % Get formatted duration string

            d = self.dur2int(self.raw_duration);
            d = self.dur2tbl(d);

        end

        function d = return_duration(self)

            d = self.raw_duration;

        end

        function i = get.raw_intervals(self)

            if isempty(self.tN), i = []; return; end
            i = diff([self.t0, self.tN, self.tF]);

        end

        function i = get.intervals(self)

            if isempty(self.tN), i = self.raw_duration; return; end
            i = self.dur2int(self.raw_intervals);
            i = self.dur2tbl(i);
            i = addvars(i, (1:self.n_lapse)', 'Before', 1, 'NewVariableNames', "rank"); % Insert before the first variable (index 1)

        end

        function d = return_intervals(self)

            d = self.raw_intervals;

        end

        function i = get.raw_lags(self)

            if isempty(self.tN), i = []; return; end
            i = [self.tN, self.tF] - self.t0;

        end

        function i = get.lags(self)

            if isempty(self.tN), i = self.raw_duration; return; end
            i = self.dur2int(self.raw_lags);
            i = self.dur2tbl(i);
            i = addvars(i, (1:self.n_lapse)', 'Before', 1, 'NewVariableNames', "rank"); % Insert before the first variable (index 1)

        end

        function d = return_lags(self)

            d = self.raw_lags;

        end

        function n = get.n_lapse(self)

            if isempty(self.tN), n = 0; return; end
            n = length(self.tN) + 1;

        end
        function varargout = report(self, type)

            if nargin == 1

                if ~isempty(self.tN)
                    type = 'intervals';
                else, type = 'duration';
                end

            elseif ~any(strcmp(type, ["intervals", "lags", "duration"]))
                error("Unrecognized type '%s' for report().", type)

            end

            dur = self.(type);
            types = string(dur.Properties.VariableNames);
            durs = table2array(dur);

            switch type

                case "intervals"

                    s_init = "Elapsed times between intervals are:\n\n";

                case "lags"

                    s_init = "Elapsed times relative to the timer onset are:\n\n";

                case "duration"
                    s_init = 'Elapsed time is ';
            end

            if strcmp(type, "duration")

                s = sprintf("%s %s.\n",s_init, join(arrayfun(@(x,y) sprintf("%d %s",x,y), durs, types),', '));

            else

                s = sprintf("%s%s",s_init, compose(join(compose("\t%d. %d %%s, %d %%s\n", durs(:,1), durs(:,2:end)),''),repmat(types(2:end),1, self.n_lapse)));

            end

            fprintf(s);
            if nargout == 1, varargout{1} = s; end

        end


    end

    methods(Access = protected)

        function t = dur2tbl(self, d)

            isIncld = any(d,1);

            t = array2table(d(:, isIncld), 'VariableNames', self.duration_types(isIncld));

        end

    end

    methods (Access = protected, Static)

        function varargout = dur2int(dur)

            funcs = gen.Timer.duration_types;
            durs = zeros(length(dur), length(funcs));
            for iFunc = 1:length(funcs)

                for iDur = 1:size(durs,1)

                    [durN, dur_tmp] = gen.Timer.count_dur_by_type(dur(iDur), funcs(iFunc));
                    dur(iDur) = dur_tmp;
                    durs(iDur, iFunc) = durN;

                end

            end

            varargout{1} = durs;
            if nargout > 1, varargout{2} = dur; end

        end

        function [subdur, rem] = count_dur_by_type(dur, func)

            subdur = floor(feval(func,dur)); % Extract duration (e.g. seconds)
            rem = dur - feval(func,subdur); % Subtract extracted time

        end

    end

end