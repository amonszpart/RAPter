function h = function_handle(fcn)
% FUNCTION_HANDLE        Construct function handle from string or
%                        cell array of strings
%
% FUNCTION_HANDLE constructs a <a href="matlab:help([matlabroot '/toolbox/matlab/datatypes/function_handle'])">function_handle</a> for built-in or 
% user-defined functions.
%
% When constructing a <a href="matlab:help([matlabroot
% '/toolbox/matlab/datatypes/function_handle'])">function_handle</a> directly with <a href="matlab:help('@')">@</a> or <a href="matlab:help('str2func')">str2func</a>, 
% it is only possible to succesfully evaluate the resulting handle when
% the function the handle refers to was on the MATLAB search at the time 
% the handle was created. To create valid handles to arbitrarty 
% functions possibly not on the search path, use the FUNCTION_HANDLE 
% constructor. 
%
% F = FUNCTION_HANDLE(fcn) for string or cell string [fcn] creates a
% function_handle or cell array of function handles for each of the 
% handles or strings in [fcn]. Any string in [fcn] may include the 
% path to the function of interest. An error is issued in case the 
% path information is incorrect, of the specified file is not an 
% executable MATLAB function.
%
% EXAMPLES: 
%
%    >> F = function_handle('./path/to/function/myFcn.m')
%    F = 
%         @myFcn
% 
%    >> A = function_handle(...
%               {@cos, '../dir/not/on/MATLAB/path/myFunction.m'})
%    A = 
%         @cos    @myFunction
% 
%    >> A{1}(pi)
%    ans = 
%       -1 
% 
%    >> functions(A{1})
%    ans = 
%         function: 'min'
%             type: 'simple'
%             file: ''
% 
%    >> functions(A{2})
%    ans = 
%         function: 'myFunction'
%             type: 'simple'
%             file: '/fullpath/dir/not/on/MATLAB/path/myFunction.m'
%
%  See also <a href="matlab:help([matlabroot '/toolbox/matlab/datatypes/function_handle'])">function_handle (built-in)</a>, str2func, functions.
   
    
% Please report bugs and inquiries to: 
%
% Name       : Rody P.S. Oldenhuis
% E-mail     : oldenhuis@gmail.com    (personal)
%              oldenhuis@luxspace.lu  (professional)
% Affiliation: LuxSpace Sàrl
% Licence    : BSD


% Changelog
%{
2014/March/19 (Rody Oldenhuis)
 - NEW: first version
%}


% If you find this work useful, please consider a donation:
% https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=6G3S5UYM7HJ3N
    
    %% Initialize
    
    % Quick exit
    if isa(fcn, 'function_handle')
        h = fcn; return; end
    
    if ~iscell(fcn)
        fcn = {fcn}; end
    
    % Extract everything that already is a function handle
    h    = cell(size(fcn));
    F    = cellfun('isclass', fcn, 'function_handle');
    h(F) = fcn(F);
    if all(F)
        return; end
    
    % Continue with the ones specified as strings
    fcn = fcn(~F);    
    if ~iscellstr(fcn)
        throwAsCaller(MException(...
            'function_handle:invalid_objects',...
            'Invalid types detected in input. Expected types are ''function_handle'' or ''char''.'));            
    end
    hF = cell(size(fcn));
    
    %% Get to work
    
    % Make sure we always end up where we started
    prevDir = pwd;
    clean__ = onCleanup(@(~)cd(prevDir));
        
    for ii = 1:numel(fcn)  
        
        % Valid inputs
        if any(exist(fcn{ii})==[2 3 5 6 8]) %#ok<EXIST>
                        
            [pth,name] = fileparts(fcn{ii});
                        
            % non-builtin
            if ~isempty(pth) 
                if exist(pth,'dir')==7
                    cd(pth);
                        hF{ii} = str2func(['@' name]);
                    cd(prevDir);
                else
                    throwAsCaller(MException(...
                        'function_handle:dir_not_found',...
                        'Directory ''%s'' not found.', pth));
                end
                
            % builtin
            elseif exist(fcn{ii},'builtin')==5
                hF{ii} = str2func(fcn{ii});
                
            % unrecognized
            else
                throwAsCaller(MException(...
                    'function_handle:fcn_not_found',...
                    ['Function ''%s'' is not on the MATLAB search path, and does not seem ',...
                    'to be a builtin.'], name));
            end
        
        % Invalid input
        else
            throwAsCaller(MException(...
                'function_handle:fcn_invalid',...
                'Function or file ''%s'' not found.', fcn{ii}));
        end
    end
    
    % Final assignment
    h(~F) = hF;
    
    % Make output more intuitive
    if numel(h)==1
        h = h{1}; end
    
end
