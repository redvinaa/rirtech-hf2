function s:Rplc()
	:s/\\mathrm{ke}/\\ke/g
	:s/\\mathrm{km}/\\km/g
	:s/\\mathrm{Ja}/\\Ja/g
	:s/\\mathrm{La}/\\La/g
	:s/\\mathrm{Ra}/\\Ra/g
	:s/\\mathrm{ks}/\\ks/g
	:s/\\mathrm{P}/P/g
	:s/\\mathrm{TI}/T_\\text{I}/g
	:s/\\mathrm{TD}/T_\\text{D}/g
	:s/\\mathrm{T0}/T_0/g
	:s/\\mathrm{tau0}/\\tau_0/g
	:s/\\mathrm{w0}/\\omega_0/g
	:s/\\mathrm{wn}/\\omega_\text{n}/g
endfunction
function s:Rplc2()
	:%s/\\omega_\\text/\\Omega_\\text/g
	:%s/u_0/U_0/g
	:%s/u_\\text{be}/U_\\text{be}/g
	:%s/i_\\text{a}/I_\\text{a}/g
	:%s/\\fn{W}_\\text{o}/\\fn{W}_\\text{tmp}/g
	:%s/\\fn{W}_\\text{x}/\\fn{W}_\\text{o}/g
	:%s/\\fn{W}_\\text{tmp}/\\fn{W}_\\text{x}/g
	:%s/\\fn{W}_\\text{f}/\\fn{W}_\\text{f}/g " TODO
	:%s/y(/\\fn{y}(/g
	:%s/I_\\text/\\fn{I}_\\text/g
	:%s/U_\\text/\\fn{U}_\\text/g
	:%s/U_0/\\fn{U}_0/g
endfunction
function s:Rplc3()
	:%s/Wo/Wtmp/g
	:%s/Wx/Wo/g
	:%s/Wtmp/Wx/g
endfunction

vnoremap <leader>v :call <SID>Rplc()<cr>
nnoremap <leader>V :call <SID>Rplc2()<cr>
nnoremap <leader>v :call <SID>Rplc3()<cr>

