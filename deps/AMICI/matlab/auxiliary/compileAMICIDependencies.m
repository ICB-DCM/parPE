function [objectsstr, includesstr] = compileAMICIDependencies(dependencyPath, objectFolder, o_suffix, COPT, DEBUG)
    %COMPILEAMICIDEPENDENCIES Compiles Sundials and SuiteSparse libraries required by AMICI
    COPT = ['CFLAGS=''$CFLAGS -std=c99'' ' COPT];
    sundials_path = fullfile(dependencyPath,'sundials');
    sundials_ver = '5.2.0';

    ssparse_path = fullfile(dependencyPath,'SuiteSparse');
    ssparse_ver = '5.4.0';

    lapack_path = fullfile(dependencyPath,'lapack-3.5.0'); % currently not used, lapack implementation still needs to be done
    lapack_ver = '3.5.0';


    version_file = fullfile(objectFolder, 'versions.txt');
    [del_sundials, del_ssparse, del_lapack] = checkVersions(version_file, sundials_ver, ssparse_ver, lapack_ver);

    % assemble objectsstr
    objectsstr = '';
    objects_sundials = getObjectsSundials(o_suffix);
    for j=1:length(objects_sundials)
        objectsstr = strcat(objectsstr,' "',fullfile(objectFolder,objects_sundials{j}),'"');
    end

    objects_ssparse = getObjectsSSparse(o_suffix);
    for j=1:length(objects_ssparse)
        objectsstr = strcat(objectsstr,' "',fullfile(objectFolder,objects_ssparse{j}),'"');
    end

    includesstr = getIncludeString(fullfile(fileparts(dependencyPath)), sundials_path, ssparse_path);

    % collect files that need to be recompiled
    sources_sundials = getSourcesSundials();
    sourcesToCompile = '';
    for j=1:length(sources_sundials)
        if(del_sundials || ~exist(fullfile(objectFolder,objects_sundials{j}), 'file'))
            sourcesToCompile = [sourcesToCompile, ' "', fullfile(sundials_path,sources_sundials{j}), '"'];
        end
    end
    sources_ssparse = getSourcesSSparse();
    for j=1:length(sources_ssparse)
        if(del_ssparse || ~exist(fullfile(objectFolder,objects_ssparse{j}), 'file'))
            sourcesToCompile = [sourcesToCompile, ' "', fullfile(ssparse_path,sources_ssparse{j}), '"'];
        end
    end

    % compile
    if(~strcmp(sourcesToCompile, ''))
        cmd = ['mex ' DEBUG ' ' COPT ' -c -outdir "' ...
            objectFolder '" ' ...
            includesstr ' ' sourcesToCompile];
        try
            eval(cmd);
        catch ME
            disp(cmd);
            rethrow(ME);
        end
    end

    % only write versions.txt if we are done compiling
    fid = fopen(version_file,'w');
    fprintf(fid,[sundials_ver '\r']);
    fprintf(fid,[ssparse_ver '\r']);
    fprintf(fid,[lapack_ver '\r']);
    fclose(fid);
end

function includesstr = getIncludeString(amici_root_path, sundials_path, ssparse_path)
    includesstr = '';
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'include', 'sundials', 'amici_matlab'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'src'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(sundials_path, 'src', 'sundials'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(amici_root_path), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(amici_root_path, 'src'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(amici_root_path, 'include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'KLU','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'AMD','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'COLAMD','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'BTF','Include'), '"');
    includesstr = strcat(includesstr,' -I"', fullfile(ssparse_path, 'SuiteSparse_config'), '"');
end


function [del_sundials, del_ssparse, del_lapack] = checkVersions(version_file, sundials_ver, ssparse_ver, lapack_ver)
    % read version number from versions.txt and decide whether object have to
    % be regenerated
    if(~exist(version_file,'file'))
        del_sundials = true;
        del_ssparse = true;
        del_lapack = true;
        fid = fopen(version_file,'w');
        fprintf(fid,[sundials_ver '\r']);
        fprintf(fid,[ssparse_ver '\r']);
        fprintf(fid,[lapack_ver '\r']);
        fclose(fid);
    else
        fid = fopen(version_file,'r');
        sundials_objver = fgetl(fid);
        ssparse_objver = fgetl(fid);
        lapack_objver = fgetl(fid);
        fclose(fid);
        del_sundials = isnewer(sundials_ver,sundials_objver);
        if(del_sundials)
            display('Newer version of Sundials! Recompiling ...')
        end
        del_ssparse = isnewer(ssparse_ver,ssparse_objver);
        if(del_ssparse)
            display('Newer version of SuiteSparse! Recompiling ...')
        end
        del_lapack = isnewer(lapack_ver,lapack_objver);
        if(del_lapack)
            display('Newer version of Lapack! Recompiling ...')
        end
    end
end

function sources_sundials = getSourcesSundials()
    sources_sundials = {
        fullfile('src', 'sunmatrix', 'dense',   'sunmatrix_dense.c');
        fullfile('src', 'sunmatrix', 'sparse',  'sunmatrix_sparse.c');
        fullfile('src', 'sunmatrix', 'band',    'sunmatrix_band.c');
        fullfile('src', 'sunlinsol', 'spgmr',   'sunlinsol_spgmr.c');
        fullfile('src', 'sunlinsol', 'sptfqmr', 'sunlinsol_sptfqmr.c');
        fullfile('src', 'sunlinsol', 'klu',     'sunlinsol_klu.c');
        fullfile('src', 'sunlinsol', 'dense',   'sunlinsol_dense.c');
        fullfile('src', 'sunlinsol', 'spfgmr',  'sunlinsol_spfgmr.c');
        fullfile('src', 'sunlinsol', 'pcg',     'sunlinsol_pcg.c');
        fullfile('src', 'sunlinsol', 'spbcgs',  'sunlinsol_spbcgs.c');
        fullfile('src', 'sunlinsol', 'band',    'sunlinsol_band.c');
        fullfile('src', 'idas', 'idaa.c');
        fullfile('src', 'idas', 'idas_ic.c');
        fullfile('src', 'idas', 'idas_nls_stg.c');
        fullfile('src', 'idas', 'idas.c');
        fullfile('src', 'idas', 'idas_bbdpre.c');
        fullfile('src', 'idas', 'idas_nls.c');
        fullfile('src', 'idas', 'idas_ls.c');
        fullfile('src', 'idas', 'idas_io.c');
        fullfile('src', 'idas', 'idas_nls_sim.c');
        fullfile('src', 'idas', 'idaa_io.c');
        fullfile('src', 'sundials', 'sundials_math.c');
        fullfile('src', 'sundials', 'sundials_matrix.c');
        fullfile('src', 'sundials', 'sundials_direct.c');
        fullfile('src', 'sundials', 'sundials_nvector_senswrapper.c');
        fullfile('src', 'sundials', 'sundials_dense.c');
        fullfile('src', 'sundials', 'sundials_nvector.c');
        fullfile('src', 'sundials', 'sundials_version.c');
        fullfile('src', 'sundials', 'sundials_iterative.c');
        fullfile('src', 'sundials', 'sundials_nonlinearsolver.c');
        fullfile('src', 'sundials', 'sundials_linearsolver.c');
        fullfile('src', 'sundials', 'sundials_band.c');
        fullfile('src', 'sundials', 'sundials_context.c');
        fullfile('src', 'sundials', 'sundials_logger.c');
        fullfile('src', 'sundials', 'sundials_errors.c');
        fullfile('src', 'sundials', 'sundials_hashmap.c');
        fullfile('src', 'sunnonlinsol', 'newton', 'sunnonlinsol_newton.c');
        fullfile('src', 'sunnonlinsol', 'fixedpoint', ...
                     'sunnonlinsol_fixedpoint.c');
        fullfile('src', 'nvector', 'serial', 'nvector_serial.c');
        fullfile('src', 'cvodes', 'cvodes_nls_stg.c');
        fullfile('src', 'cvodes', 'cvodes_ls.c');
        fullfile('src', 'cvodes', 'cvodes_nls_stg1.c');
        fullfile('src', 'cvodes', 'cvodes_bbdpre.c');
        fullfile('src', 'cvodes', 'cvodes.c');
        fullfile('src', 'cvodes', 'cvodes_bandpre.c');
        fullfile('src', 'cvodes', 'cvodea.c');
        fullfile('src', 'cvodes', 'cvodes_nls_sim.c');
        fullfile('src', 'cvodes', 'cvodea_io.c');
        fullfile('src', 'cvodes', 'cvodes_nls.c');
        fullfile('src', 'cvodes', 'cvodes_diag.c');
        fullfile('src', 'cvodes', 'cvodes_io.c');
        fullfile('src', 'cvodes', 'cvodes_proj.c');
        };
end



function sources_ssparse = getSourcesSSparse()
    sources_ssparse = {
        fullfile('KLU','Source','klu_l_analyze_given.c');
        fullfile('KLU','Source','klu_l_analyze.c');
        fullfile('KLU','Source','klu_l_defaults.c');
        fullfile('KLU','Source','klu_l_diagnostics.c');
        fullfile('KLU','Source','klu_l_dump.c');
        fullfile('KLU','Source','klu_l_extract.c');
        fullfile('KLU','Source','klu_l_factor.c');
        fullfile('KLU','Source','klu_l_free_numeric.c');
        fullfile('KLU','Source','klu_l_free_symbolic.c');
        fullfile('KLU','Source','klu_l_kernel.c');
        fullfile('KLU','Source','klu_l_memory.c');
        fullfile('KLU','Source','klu_l_refactor.c');
        fullfile('KLU','Source','klu_l_scale.c');
        fullfile('KLU','Source','klu_l_sort.c');
        fullfile('KLU','Source','klu_l_solve.c');
        fullfile('KLU','Source','klu_l_tsolve.c');
        fullfile('KLU','Source','klu_l.c');
        fullfile('AMD','Source','amd_l1.c');
        fullfile('AMD','Source','amd_l2.c');
        fullfile('AMD','Source','amd_l_aat.c');
        fullfile('AMD','Source','amd_l_control.c');
        fullfile('AMD','Source','amd_l_defaults.c');
        fullfile('AMD','Source','amd_l_dump.c');
        fullfile('AMD','Source','amd_l_info.c');
        fullfile('AMD','Source','amd_l_order.c');
        fullfile('AMD','Source','amd_l_postorder.c');
        fullfile('AMD','Source','amd_l_post_tree.c');
        fullfile('AMD','Source','amd_l_preprocess.c');
        fullfile('AMD','Source','amd_l_valid.c');
        fullfile('COLAMD','Source','colamd_l.c');
        fullfile('BTF','Source','btf_l_maxtrans.c');
        fullfile('BTF','Source','btf_l_order.c');
        fullfile('BTF','Source','btf_l_strongcomp.c');
        fullfile('SuiteSparse_config','SuiteSparse_config.c');
        };
end

function objects_ssparse = getObjectsSSparse(o_suffix)

    objects_ssparse = {
        'klu_l_analyze_given.o';
        'klu_l_analyze.o';
        'klu_l_defaults.o';
        'klu_l_diagnostics.o';
        'klu_l_dump.o';
        'klu_l_extract.o';
        'klu_l_factor.o';
        'klu_l_free_numeric.o';
        'klu_l_free_symbolic.o';
        'klu_l_kernel.o';
        'klu_l_memory.o';
        'klu_l_refactor.o';
        'klu_l_scale.o';
        'klu_l_sort.o';
        'klu_l_solve.o';
        'klu_l_tsolve.o';
        'klu_l.o';
        'amd_l1.o';
        'amd_l2.o';
        'amd_l_aat.o';
        'amd_l_control.o';
        'amd_l_defaults.o';
        'amd_l_dump.o';
        'amd_l_info.o';
        'amd_l_order.o';
        'amd_l_post_tree.o';
        'amd_l_postorder.o';
        'amd_l_preprocess.o';
        'amd_l_valid.o';
        'colamd_l.o';
        'btf_l_maxtrans.o';
        'btf_l_order.o';
        'btf_l_strongcomp.o';
        'SuiteSparse_config.o';
        };

    if(~strcmp(o_suffix, '.o'))
        objects_ssparse = strrep(objects_ssparse, '.o', o_suffix);
    end
end

function objects_sundials = getObjectsSundials(o_suffix)
    objects_sundials = {
        'sunmatrix_dense.o';
        'sunlinsol_spgmr.o';
        'sunlinsol_sptfqmr.o';
        'sunlinsol_klu.o';
        'idaa.o';
        'idas_ic.o';
        'idas_nls_stg.o';
        'idas.o';
        'idas_bbdpre.o';
        'idas_nls.o';
        'idas_ls.o';
        'idas_io.o';
        'idas_nls_sim.o';
        'idaa_io.o';
        'sundials_context.o';
        'sundials_logger.o';
        'sundials_errors.o';
        'sundials_hashmap.o';
        'sundials_math.o';
        'sundials_matrix.o';
        'sundials_direct.o';
        'sundials_nvector_senswrapper.o';
        'sundials_dense.o';
        'sundials_nvector.o';
        'sundials_version.o';
        'sundials_iterative.o';
        'sundials_nonlinearsolver.o';
        'sundials_linearsolver.o';
        'sundials_band.o';
        'sunmatrix_band.o';
        'sunmatrix_sparse.o';
        'sunnonlinsol_newton.o';
        'sunnonlinsol_fixedpoint.o';
        'nvector_serial.o';
        'sunlinsol_pcg.o';
        'sunlinsol_dense.o';
        'sunlinsol_spbcgs.o';
        'sunlinsol_band.o';
        'sunlinsol_spfgmr.o';
        'cvodes_nls_stg.o';
        'cvodes_ls.o';
        'cvodes_nls_stg1.o';
        'cvodes_bbdpre.o';
        'cvodes.o';
        'cvodes_bandpre.o';
        'cvodea.o';
        'cvodes_nls_sim.o';
        'cvodea_io.o';
        'cvodes_nls.o';
        'cvodes_diag.o';
        'cvodes_io.o';
        'cvodes_proj.o';
        };

    if(~strcmp(o_suffix, '.o'))
        objects_sundials = strrep(objects_sundials, '.o', o_suffix);
    end
end

function result = isnewer(ver1str,ver2str)
    % isnewer checks whether the version indicated in ver1str is newer than
    % the on in ver2str
    %
    % Parameters:
    %  ver1str: version string 1, this should be a string in the format %d.%d.%d @type string
    %  ver2str: version string 1, this should be a string in the format %d.%d.%d @type string
    %
    % Return values:
    %  result: flag indicating whether ver1 ist newer than ver2 @type boolean

    ver1Parts = getParts(ver1str);
    ver2Parts = getParts(ver2str);
    if ver2Parts(1) ~= ver1Parts(1)     % major version
        result = ver2Parts(1) < ver1Parts(1);
    elseif ver2Parts(2) ~= ver1Parts(2) % minor version
        result = ver2Parts(2) < ver1Parts(2);
    elseif ver2Parts(3) ~= ver1Parts(3) % revision version
        result = ver2Parts(3) < ver1Parts(3);
    else
        result = ver2Parts(4) < ver1Parts(4);
    end
end

function parts = getParts(V)
    % getParts takes an input version string and returns an array
    % containing the major minor and subminor version number
    %
    % Parameters:
    %  V: version string, this should be a string in the format %d.%d.%d @type string
    %
    % Return values:
    %  parts: array containing the version numbers @type double

    parts = sscanf(V, '%d.%d.%d.%d')';
    if length(parts) < 3
        parts(3) = 0; % zero-fills to 3 elements
    end
    if length(parts) < 4
        parts(4) = 0; % zero-fills to 3 elements
    end
end
