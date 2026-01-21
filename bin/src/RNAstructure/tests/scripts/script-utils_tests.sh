# Unit tests for script-utils.sh
SCRIPT_DIR=$(dirname "${BASH_SOURCE[0]}")
source "${SCRIPT_DIR}"/script-utils.sh

function test_all() {
    test_joinArgs
    test_cleanpath
    test_relpath
}

function test_cleanpath() {
    local tests
    tests=(
        '..'='..'
        '../'='../'
        '.'='.'
        './'='./'
        'dir'='dir'
        '/dir'='/dir'
        'dir/'='dir/'
        '/dir/'='/dir/'

        '../..'='../..'
        '../../'='../../'
        './..'='..'
        './../'='../'
        './dir'='dir'
        './dir/'='dir/'
        './dir/abc'='dir/abc'
        './../dir'='../dir'
        './../dir/'='../dir/'
        './../dir/abc'='../dir/abc'

        './dir/..'='.'
        './dir/../'='./'
        './dir/../abc'='abc'
        './../dir/..'='..'
        './../dir/../'='../'
        './../dir/../..'='../..'
        './../dir/../../'='../../'

        'dir/..'='.'
        'dir/../'='./'
        'dir/../abc'='abc'
        '../dir/..'='..'
        '../dir/../'='../'
        '../dir/../..'='../..'
        '../dir/../../'='../../'

        'dir//abc//123'='dir/abc/123'
        'dir//abc//123//'='dir/abc/123/'
        './dir//abc//123//'='dir/abc/123/'
        './dir/./abc/./123/./'='dir/abc/123/'
        './dir/./abc/./123/.'='dir/abc/123'
        
        'dir/../abc//123'='abc/123'
        '/dir/../abc//123'='/abc/123'
        'dir/abc/../123'='dir/123'
        '/dir/./abc/./.././123'='/dir/123'
        '/dir/./abc/././123/..'='/dir/abc'
        '/dir/./abc/././123/../'='/dir/abc/'
        '/dir/./abc/././123/./.././'='/dir/abc/'
        '/dir/./abc/././123/./../.'='/dir/abc'
        
        '/dir/abc/../../123'='/123'
        'dir/abc/../../123'='123'
        '/dir/./abc/./.././../123'='/123'
        'dir/./abc/./.././../123'='123'

        '/dir/abc/../../../123'='/../123'
        'dir/abc/../../../123'='../123'
        '/dir/./abc/./.././.././.././123'='/../123'
        'dir/./abc/./.././.././../123'='../123'

        '/dir/abc/../../123/..'='/'
        'dir/abc/../../123/..'='.'
        '/dir/./abc/./.././../123/..'='/'
        'dir/./abc/./.././../123/..'='.'
        '/dir/abc/../../123/../'='/'
        'dir/abc/../../123/../'='./'
        '/dir/./abc/./.././../123/../'='/'
        'dir/./abc/./.././../123/../'='./'

        '/dir/abc/../../123/../..'='/..'
        'dir/abc/../../123/../..'='..'
        '/dir/./abc/./.././../123/../..'='/..'
        'dir/./abc/./.././../123/../..'='..'
        '/dir/abc/../../123/../../'='/../'
        'dir/abc/../../123/../../'='../'
        '/dir/./abc/./.././../123/../../'='/../'
        'dir/./abc/./.././../123/../../'='../'

        '/dir/abc/../../../123/..'='/..'
        'dir/abc/../../../123/..'='..'
        '/dir/./abc/./.././../../123/..'='/..'
        'dir/./abc/./.././../../123/..'='..'
        '/dir/abc/../../../123/../'='/../'
        'dir/abc/../../../123/../'='../'
        '/dir/./abc/./.././../../123/../'='/../'
        'dir/./abc/./.././../../123/../'='../'
    )
    local t=$(($(date +%s%N)/1000000))
    local test arg val
    for test in "${tests[@]}"; do
        arg="${test%=*}"
        exp="${test#*=}"
        cleanpath "$arg" --ret:val
        if [[ "$val" != "$exp" ]]; then
           echo >&2 "FAILED cleanpath test for '$arg' -- expected: '$exp'  actual: '$val'"
           return 1
        fi
    done
    echo "cleanpath tests PASSED! " $(($(date +%s%N)/1000000-t)) "ms"
}

function test_relpath() {
    local tests
    tests=(
        # 'subject','reference'='result'
        'AAA/BBB','AAA'='BBB'
    )
    local t=$(($(date +%s%N)/1000000))
    local test args val
    for test in "${tests[@]}"; do
        IFS=',' read -ra args <<< "${test%=*}"  # arguments are separated from each other by a comma and from the expected result by an equals 
        exp="${test#*=}"
        relpath "${args[@]}" --ret:val
        if [[ "$val" != "$exp" ]]; then
            argstr=$(joinArgs "','" "${args[@]}")
           echo >&2 "FAILED relpath test for ['$argstr'] -- expected: '$exp'  actual: '$val'"
           return 1
        fi
    done
    echo "relpath Tests PASSED! " $(($(date +%s%N)/1000000-t)) "ms"
}

function test_joinArgs() {
    local tests
    tests=(
        # 'subject','reference'='result'
        '+,A'='A'
        '+,A,B'='A+B'
        '+,A,B,C'='A+B+C'
        '--,A,B,C'='A--B--C'
    )
    local t=$(($(date +%s%N)/1000000))
    local test args val
    for test in "${tests[@]}"; do
        IFS=',' read -ra args <<< "${test%=*}"  # arguments are separated from each other by a comma and from the expected result by an equals 
        exp="${test#*=}"
        joinArgs --ret:__val "${args[@]}"
        if [[ "$val" != "$exp" ]]; then
            echo >&2 "FAILED joinArgs test for ['${test%=*}'] -- expected: '$exp'  actual: '$val'"
            local arg i=0
            for arg in "${args[@]}"; do
              echo >&2 "    $i = ${args[i]}"
              ((++i))
            done
           return 1
        fi
    done
    echo "joinArgs Tests PASSED! " $(($(date +%s%N)/1000000-t)) "ms"
}
