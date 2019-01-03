function  ScreenComment(Comment1,Comment2)
%ScreenComment Display desired level of screen comments.
%   ScreenComment(COMMENT1,COMMENT2) will, 
%   depending on value of global CMNT_LEVEL display:
%   0: Display nothing
%   1: Display concise comments stored in COMMENT1 string
%   2: Extensive comments stored in COMMENT2 string

global cmnt_level
if cmnt_level == 1
    disp(Comment1)
elseif cmnt_level == 2
    disp(Comment2)
else
    % Display nothing
end

end

