function [e_term,s_term] = sum_fun(sigma_100,f,an,n)

fn=(125*sigma_100)^(0.8312)*10^(n-1);

e_term=an./(1+(f./fn).^2);
s_term=(an.*f./fn)./(1+(f./fn).^2);