#!/usr/bin/node

let s = 'AAGTCGUATCGTA';
if (process.argv.length < 3) {
	//console.error('Error: missing sequence argument.');
	//process.exit(1);
} else
	s = process.argv[2];

console.log('Sequence: \'%s\'', s);
let m = new Map([
	['A','U'],
	['T','A'],
	['U','A'],
	['C','G'],
	['G','C'],
]);

//for (let k of s)
//	console.log(k,'=', m.get(k));

console.log(s.split('').reduce((agg, cur)=>agg+=m.get(cur), '').split("").reverse().join(""));