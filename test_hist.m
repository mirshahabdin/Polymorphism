function v = test_hist(t)

v = zeros(1, 3);

if (t == 0)
    v(1) = 0.0;
    v(2) = 1.0;
    v(3) = 1.0;
end

% if (t > -28.0)
%     h = dde23('test_hist_coev', 28.0, 'test_hist2', [0, 28+t]);
%     v(1) = h.y(1, end);
%     v(3) = h.y(2, end);
% end
% if (t == -28)
%     v(1) = 15.0;
%     v(3) = 1.0;
% end
% if (t == 0)  
%     v(2) = 15.0 * v(3);
% end


end