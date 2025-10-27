'use client'

import React, { useState, useMemo } from 'react';
import { Upload, Play, Download, TrendingUp, AlertCircle } from 'lucide-react';
import { LineChart, Line, XAxis, YAxis, CartesianGrid, Tooltip, Legend, ResponsiveContainer } from 'recharts';

const RRNAAnalysis = () => {
  const [uploadedFiles, setUploadedFiles] = useState<any[]>([]);
  const [strainType, setStrainType] = useState('normal');
  const [simResults, setSimResults] = useState<any[]>([]);
  const [observedData, setObservedData] = useState(null);
  const [isSimulating, setIsSimulating] = useState(false);
  const [numSimulations, setNumSimulations] = useState(100);
  const [initialDivergence, setInitialDivergence] = useState(0.002);

  const params = {
    normal: {
      mutationRate: 3e-6 / 30000,
      label: 'Normal (μ = 1×10⁻¹⁰ per bp)',
      geneConvRate: 8.6e-6
    },
    mutator: {
      mutationRate: (3e-6 / 30000) * 50,
      label: 'Mutator (μ = 5×10⁻⁹ per bp)',
      geneConvRate: 8.6e-6
    }
  };

  const OPERON_LENGTHS = {
    rrsB: 1541,
    rrlB: 2904
  };

  const handleFileUpload = async (e: React.ChangeEvent<HTMLInputElement>) => {
    const files = Array.from(e.target.files || []);
    
    for (const file of files) {
      const text = await file.text();
      const lines = text.trim().split('\n');
      
      const filename = file.name;
      const match = filename.match(/(Ara[+-]\d+)_(\d+)gen_(\w+)/i);
      
      if (!match) {
        console.warn(`Could not parse filename: ${filename}`);
        continue;
      }
      
      const strain = match[1];
      const generation = parseInt(match[2]);
      const fileId = match[3];
      
      const data = lines.slice(1).map(line => {
        const [seq_id, position, ref_base, new_base, new_cov, ref_cov, total_cov, frequency] = line.split(',');
        return {
          seq_id: seq_id?.trim(),
          position: parseInt(position),
          ref_base: ref_base?.trim(),
          new_base: new_base?.trim(),
          new_cov: parseInt(new_cov),
          ref_cov: parseInt(ref_cov),
          total_cov: parseInt(total_cov),
          frequency: parseFloat(frequency)
        };
      }).filter(d => d.seq_id && !isNaN(d.position));
      
      setUploadedFiles(prev => [...prev, {
        filename: file.name,
        strain,
        generation,
        fileId,
        data
      }]);
    }
  };

  const processObservedData = useMemo(() => {
    if (uploadedFiles.length === 0) return null;
  
    const strainGenData: any = {};
    
    uploadedFiles.forEach(file => {
      const key = `${file.strain}_${file.generation}`;
      if (!strainGenData[key]) {
        strainGenData[key] = {
          strain: file.strain,
          generation: file.generation,
          files: [],
          allData: []
        };
      }
      strainGenData[key].files.push(file.fileId);
      strainGenData[key].allData.push(...file.data);
    });
  
    const strainData: any = {};
    
    Object.values(strainGenData).forEach((entry: any) => {
      if (!strainData[entry.strain]) {
        strainData[entry.strain] = [];
      }
      
      const uniquePolymorphisms = new Map();
      entry.allData.forEach((d: any) => {
        const key = `${d.seq_id}_${d.position}`;
        if (!uniquePolymorphisms.has(key) || uniquePolymorphisms.get(key).frequency < d.frequency) {
          uniquePolymorphisms.set(key, d);
        }
      });
      
      const polymorphisms = Array.from(uniquePolymorphisms.values());
  
      const OPERONLENGTHS = { rrsB: 1541, rrlB: 2904 };
  
      const rrsmutations = polymorphisms.filter(d => d.seq_id === "rrsB");
      const rrlmutations = polymorphisms.filter(d => d.seq_id === "rrlB");
  
      // Divergence: frequency-weighted, normalized by all possible sites!
      const sumFreqRrs = rrsmutations.reduce((sum, d) => sum + d.frequency, 0);
      const sumFreqRrl = rrlmutations.reduce((sum, d) => sum + d.frequency, 0);
  
      const rrs_divergence = sumFreqRrs / OPERONLENGTHS.rrsB;
      const rrl_divergence = sumFreqRrl / OPERONLENGTHS.rrlB;
  
      // Weighted average over all sites of both operons
      const avg_divergence = (sumFreqRrs + sumFreqRrl) / (OPERONLENGTHS.rrsB + OPERONLENGTHS.rrlB);
      
      strainData[entry.strain].push({
        generation: entry.generation,
        avg_divergence: avg_divergence,
        rrs_divergence: rrs_divergence,
        rrl_divergence: rrl_divergence,
        total_polymorphisms: polymorphisms.length,
        rrs_count: rrsmutations.length,
        rrl_count: rrlmutations.length,
        numFiles: entry.files.length
      });
    });
  
    Object.keys(strainData).forEach(strain => {
      strainData[strain].sort((a: any, b: any) => a.generation - b.generation);
    });
  
    return strainData;
  }, [uploadedFiles]);
  

  const runSimulation = (numRuns = 100) => {
    setIsSimulating(true);
    
    setTimeout(() => {
      const results = [];
      const currentParams = params[strainType as keyof typeof params];
      const operons = ['rrnA', 'rrnB', 'rrnC', 'rrnD', 'rrnE', 'rrnG', 'rrnH'];
      
      for (let run = 0; run < numRuns; run++) {
        let rrsSeqs: any = {}, rrlSeqs: any = {};
        operons.forEach(op => {
          rrsSeqs[op] = new Array(OPERON_LENGTHS.rrsB).fill(0);
          rrlSeqs[op] = new Array(OPERON_LENGTHS.rrlB).fill(0);
          
          if (op !== 'rrnB') {
            const totalSites = OPERON_LENGTHS.rrsB + OPERON_LENGTHS.rrlB;
            const numDiffSites = Math.floor(totalSites * initialDivergence);
            
            const diffSitesInRrs = Math.floor(numDiffSites * (OPERON_LENGTHS.rrsB / totalSites));
            const diffSitesInRrl = numDiffSites - diffSitesInRrs;
            
            for (let i = 0; i < diffSitesInRrs; i++) {
              const pos = Math.floor(Math.random() * OPERON_LENGTHS.rrsB);
              rrsSeqs[op][pos] = 1 + Math.floor(Math.random() * 3);
            }
            
            for (let i = 0; i < diffSitesInRrl; i++) {
              const pos = Math.floor(Math.random() * OPERON_LENGTHS.rrlB);
              rrlSeqs[op][pos] = 1 + Math.floor(Math.random() * 3);
            }
          }
        });
        
        const trajectory = [];
        
        for (let gen = 0; gen <= 50000; gen++) {
          operons.forEach(operon => {
            if (Math.random() < (currentParams.mutationRate * OPERON_LENGTHS.rrsB)) {
              const pos = Math.floor(Math.random() * OPERON_LENGTHS.rrsB);
              rrsSeqs[operon][pos] = (rrsSeqs[operon][pos] + 1) % 4;
            }
            
            if (Math.random() < (currentParams.mutationRate * OPERON_LENGTHS.rrlB)) {
              const pos = Math.floor(Math.random() * OPERON_LENGTHS.rrlB);
              rrlSeqs[operon][pos] = (rrlSeqs[operon][pos] + 1) % 4;
            }
          });
          
          operons.forEach(recipient => {
            operons.forEach(donor => {
              if (donor !== recipient && Math.random() < currentParams.geneConvRate) {
                const tract_len = 50 + Math.floor(Math.random() * 150);
                const start = Math.floor(Math.random() * Math.max(1, OPERON_LENGTHS.rrsB - tract_len));
                for (let i = start; i < Math.min(start + tract_len, OPERON_LENGTHS.rrsB); i++) {
                  rrsSeqs[recipient][i] = rrsSeqs[donor][i];
                }
                
                const start2 = Math.floor(Math.random() * Math.max(1, OPERON_LENGTHS.rrlB - tract_len));
                for (let i = start2; i < Math.min(start2 + tract_len, OPERON_LENGTHS.rrlB); i++) {
                  rrlSeqs[recipient][i] = rrlSeqs[donor][i];
                }
              }
            });
          });
          
          if (gen % 2500 === 0) {
            const refRrs = rrsSeqs['rrnB'];
            const refRrl = rrlSeqs['rrnB'];
            
            let totalDivergence = 0;
            let numOperons = 0;
            
            operons.forEach(op => {
              if (op !== 'rrnB') {
                let diffSites = 0;
                
                for (let i = 0; i < OPERON_LENGTHS.rrsB; i++) {
                  if (rrsSeqs[op][i] !== refRrs[i]) {
                    diffSites++;
                  }
                }
                
                for (let i = 0; i < OPERON_LENGTHS.rrlB; i++) {
                  if (rrlSeqs[op][i] !== refRrl[i]) {
                    diffSites++;
                  }
                }
                
                const operonDiv = diffSites / (OPERON_LENGTHS.rrsB + OPERON_LENGTHS.rrlB);
                totalDivergence += operonDiv;
                numOperons++;
              }
            });
            
            const avgDiv = numOperons > 0 ? totalDivergence / numOperons : 0;
            
            trajectory.push({
              generation: gen,
              total_divergence: avgDiv
            });
          }
        }
        
        results.push(trajectory);
      }
      
      setSimResults(results);
      setIsSimulating(false);
    }, 100);
  };

  type SimResultItem = {
    generation: number;
    total_divergence: number;
  };
  
  const summaryStats = useMemo(() => {
    if (simResults.length === 0) return null;
  
    // Explicitly type d as SimResultItem
    const generations = simResults[0].map((d: SimResultItem) => d.generation);
  
    const summary = generations.map((gen: number) => {
      // Explicitly type run as array of SimResultItem
      const values = simResults
        .map((run: SimResultItem[]) =>
          run.find((d: SimResultItem) => d.generation === gen)?.total_divergence || 0
        )
        // Type guard to tell TS these are numbers
        .filter((v): v is number => v !== undefined && !isNaN(v));
  
      // Sort values for median and quartiles
      values.sort((a, b) => a - b);
  
      const mean = values.reduce((a, b) => a + b, 0) / values.length;
      const median = values[Math.floor(values.length / 2)];
      const q25 = values[Math.floor(values.length * 0.25)];
      const q75 = values[Math.floor(values.length * 0.75)];
  
      return {
        generation: gen,
        mean,
        median,
        q25,
        q75,
        min: values[0],
        max: values[values.length - 1],
      };
    });
  
    return summary;
  }, [simResults]);
  

  return (
    <div className="p-6 max-w-7xl mx-auto bg-gray-50 min-h-screen">
      <div className="bg-white rounded-lg shadow-lg p-6 mb-6">
        <h1 className="text-3xl font-bold text-gray-800 mb-2">
          rRNA Operon Evolution: Neutral Model Analysis
        </h1>
        <p className="text-gray-600">
          Compare observed rRNA operon divergence (all operons vs rrnB) against neutral expectations
        </p>
      </div>

      <div className="bg-white rounded-lg shadow p-6 mb-6">
        <h2 className="text-xl font-semibold mb-4 flex items-center gap-2">
          <Upload size={20} />
          Upload Data Files
        </h2>
        <div className="space-y-4">
          <div>
            <input
              type="file"
              multiple
              accept=".csv"
              onChange={handleFileUpload}
              className="w-full text-sm"
            />
            <p className="text-xs text-gray-500 mt-2">
              Upload CSV files named like: Ara-1_500gen_762B.csv
              <br />
              Expected format: [Strain]_[Generation]gen_[FileID].csv
            </p>
          </div>
          
          <div className="bg-gray-50 p-3 rounded">
            <p className="text-sm font-medium mb-1">Files uploaded: {uploadedFiles.length}</p>
            {uploadedFiles.length > 0 && (
              <div className="text-xs text-gray-600 max-h-32 overflow-y-auto space-y-1">
                {uploadedFiles.slice(0, 10).map((f, i) => (
                  <div key={i}>
                    {f.filename} → {f.strain} @ {f.generation}gen (ID: {f.fileId})
                  </div>
                ))}
                {uploadedFiles.length > 10 && <div>... and {uploadedFiles.length - 10} more</div>}
              </div>
            )}
            {processObservedData && (
              <div className="mt-2 pt-2 border-t">
                <p className="text-xs font-medium">Data Summary:</p>
                <p className="text-xs text-gray-600">
                  {Object.keys(processObservedData).length} strains, 
                  {' '}{Object.values(processObservedData).reduce((sum: number, arr: any) => sum + arr.length, 0)} timepoints total
                </p>
              </div>
            )}
          </div>
        </div>
      </div>

      <div className="bg-white rounded-lg shadow p-6 mb-6">
        <h2 className="text-xl font-semibold mb-4 flex items-center gap-2">
          <Play size={20} />
          Run Neutral Model Simulations
        </h2>
        
        <div className="grid grid-cols-1 md:grid-cols-3 gap-4 mb-4">
          <div>
            <label className="block text-sm font-medium mb-2">Strain Type</label>
            <select
              value={strainType}
              onChange={(e) => setStrainType(e.target.value)}
              className="w-full border rounded px-3 py-2"
            >
              <option value="normal">Normal mutation rate</option>
              <option value="mutator">Mutator (50x higher μ)</option>
            </select>
          </div>
          
          <div>
            <label className="block text-sm font-medium mb-2">Number of Simulations</label>
            <input
              type="number"
              value={numSimulations}
              onChange={(e) => setNumSimulations(parseInt(e.target.value))}
              className="w-full border rounded px-3 py-2"
              min={10}
              max={1000}
            />
          </div>

          <div>
            <label className="block text-sm font-medium mb-2">
              Initial Divergence
              <span className="text-xs text-gray-500 ml-1">(0-1)</span>
            </label>
            <input
              type="number"
              step={0.001}
              value={initialDivergence}
              onChange={(e) => setInitialDivergence(parseFloat(e.target.value))}
              className="w-full border rounded px-3 py-2"
              min={0}
              max={1}
            />
            <p className="text-xs text-gray-500 mt-1">Fraction of sites differing at t=0</p>
          </div>
        </div>

        <div className="mb-4">
          <button
            onClick={() => runSimulation(numSimulations)}
            disabled={isSimulating}
            className="w-full bg-blue-600 text-white px-6 py-2 rounded hover:bg-blue-700 disabled:bg-gray-400"
          >
            {isSimulating ? 'Simulating...' : 'Run Simulations'}
          </button>
        </div>

        <div className="bg-blue-50 p-3 rounded text-sm">
          <p className="font-medium mb-1">Current Parameters:</p>
          <ul className="text-xs space-y-1">
            <li>• Mutation rate: {params[strainType as keyof typeof params].mutationRate.toExponential(2)} per bp per generation</li>
            <li>• Gene conversion: 8.6×10⁻⁶ per operon per donor per generation</li>
            <li>• Initial divergence: {initialDivergence} (fraction of sites different at t=0)</li>
            <li>• Generations: 0 to 50,000</li>
          </ul>
        </div>
      </div>

      {summaryStats && (
        <div className="bg-white rounded-lg shadow p-6 mb-6">
          <h2 className="text-xl font-semibold mb-4 flex items-center gap-2">
            <TrendingUp size={20} />
            Simulation Results: Neutral Expectation
          </h2>
          
          <ResponsiveContainer width="100%" height={400}>
          <LineChart data={summaryStats}
            margin={{ top: 20, right: 20, bottom: 20, left: 60 }}  // increase 'left' as needed
              >
              <CartesianGrid strokeDasharray="3 3" />
              <XAxis 
                dataKey="generation" 
                label={{ value: 'Generation', position: 'insideBottom', offset: -5,}}
              />
              <YAxis label={{ value: 'Avg Divergence from rrnB', angle: -90, position: 'insideLeft', offset: -20, dy: 60 }} />
              <Tooltip />
              <Legend />
              <Line type="monotone" dataKey="mean" stroke="#2563eb" strokeWidth={2} name="Mean" dot={false} />
              <Line type="monotone" dataKey="q25" stroke="#93c5fd" strokeWidth={1.5} name="25th percentile" dot={false} />
              <Line type="monotone" dataKey="q75" stroke="#93c5fd" strokeWidth={1.5} name="75th percentile" dot={false} />
              <Line type="monotone" dataKey="median" stroke="#1e40af" strokeWidth={1} strokeDasharray="5 5" name="Median" dot={false} />
            </LineChart>
          </ResponsiveContainer>
          
          <div className="mt-4 p-4 bg-gray-50 rounded">
            <p className="text-sm font-medium mb-2">Interpretation:</p>
            <ul className="text-sm space-y-1 text-gray-700">
              <li>• Blue line = average expected divergence under neutral evolution</li>
              <li>• Light blue lines = middle 50% of simulations (25th-75th percentile)</li>
              <li>• Dashed line = median</li>
              <li>• Gene conversion balances mutation, creating equilibrium</li>
              <li>• Compare your observed data to these bands</li>
            </ul>
            <div className="mt-2 text-xs text-gray-600">
              <p>Debug: Mean at 50k gen = {summaryStats[summaryStats.length - 1]?.mean?.toFixed(4)}, 
              Q25 = {summaryStats[summaryStats.length - 1]?.q25?.toFixed(4)}, 
              Q75 = {summaryStats[summaryStats.length - 1]?.q75?.toFixed(4)}</p>
            </div>
          </div>
        </div>
      )}

      {processObservedData && Object.keys(processObservedData).length > 0 && (
        <div className="bg-white rounded-lg shadow p-6 mb-6">
          <h2 className="text-xl font-semibold mb-4">Observed Data</h2>
          
          <div className="mb-4 p-4 bg-blue-50 rounded">
            <p className="text-sm font-medium mb-2">Data Summary:</p>
            <div className="text-xs space-y-1">
              {Object.entries(processObservedData).map(([strain, data]: [string, any]) => (
                <div key={strain}>
                  <strong>{strain}:</strong> {data.length} timepoints
                  {data.length > 0 && (
                    <span className="ml-2 text-gray-600">
                      (Gen: {Math.min(...data.map((d: any) => d.generation))} - {Math.max(...data.map((d: any) => d.generation))}, 
                      Avg Div: {Math.min(...data.map((d: any) => d.avg_divergence)).toFixed(3)} - {Math.max(...data.map((d: any) => d.avg_divergence)).toFixed(3)})
                    </span>
                  )}
                </div>
              ))}
            </div>
          </div>

          <div style={{ width: '100%', height: 400 }}>
            <ResponsiveContainer>
              <LineChart margin={{ top: 20, right: 20, bottom: 60, left: 60, }}>
                <CartesianGrid strokeDasharray="3 3" />
                <XAxis 
                  dataKey="generation" 
                  type="number"
                  domain={[0, 50000]}
                  label={{ value: 'Generation', position: 'insideBottom', offset: -5 }}
                />
                <YAxis
                  // ... other props
                  domain={[
                    'auto',
                    (dataMax: number) => Math.min(dataMax * 1.05)
                  ]}
                />

                <Tooltip 
                  content={({ active, payload }) => {
                    if (active && payload && payload.length) {
                      const data = payload[0].payload;
                      return (
                        <div className="bg-white p-2 border rounded shadow text-xs">
                          <p className="font-medium">{payload[0].name}</p>
                          <p>Generation: {data.generation}</p>
                          <p>Avg Divergence: {data.avg_divergence?.toFixed(4)}</p>
                          <p>16S: {data.rrs_divergence?.toFixed(4)} ({data.rrs_count} sites)</p>
                          <p>23S: {data.rrl_divergence?.toFixed(4)} ({data.rrl_count} sites)</p>
                          <p className="text-gray-500 mt-1">Total polymorphisms: {data.total_polymorphisms}</p>
                        </div>
                      );
                    }
                    return null;
                  }}
                />
                <Legend />
                {Object.entries(processObservedData).map(([strain, data]: [string, any], i) => {
                  const isMutator = ['Ara-2', 'Ara-4', 'Ara+3', 'Ara-1', 'Ara+6', 'Ara-3',].includes(strain);
                  
                  return (
                    <Line
                      key={strain}
                      name={strain}
                      data={data}
                      dataKey="avg_divergence"
                      stroke={isMutator ? '#ef4444' : `hsl(${(i * 137) % 360}, 70%, 50%)`}
                      strokeWidth={isMutator ? 2 : 1.5}
                      dot={{ r: 4 }}
                      type="monotone"
                      connectNulls={false}
                    />
                  );
                })}
              </LineChart>
            </ResponsiveContainer>
          </div>
          
          <div className="mt-4 grid grid-cols-2 md:grid-cols-4 gap-4">
            {Object.entries(processObservedData).map(([strain, data]: [string, any]) => {
              const isMutator = ['Ara-1', 'Ara+6', 'Ara-2', 'Ara-3', 'Ara-4', 'Ara+3', ].includes(strain);
              return (
                <div key={strain} className={`p-3 rounded ${isMutator ? 'bg-red-50 border border-red-200' : 'bg-gray-50'}`}>
                  <p className="font-medium text-sm">
                    {strain}
                    {isMutator && <span className="ml-1 text-xs text-red-600">(mutator)</span>}
                  </p>
                  <p className="text-xs text-gray-600">Timepoints: {data.length}</p>
                  <p className="text-xs text-gray-600">
                    Range: {data[0]?.generation || 0} - {data[data.length - 1]?.generation || 0} gen
                  </p>
                </div>
              );
            })}
          </div>
        </div>
      )}

      {summaryStats && processObservedData && (
        <div className="bg-white rounded-lg shadow p-6 mb-6">
          <h2 className="text-xl font-semibold mb-4 flex items-center gap-2">
            <AlertCircle size={20} />
            Statistical Comparison
          </h2>
          
          <div className="space-y-4">
            <div className="p-4 bg-yellow-50 rounded">
              <p className="text-sm font-medium mb-2">Next Steps for Formal Testing:</p>
              <ol className="text-sm space-y-1 list-decimal list-inside">
                <li>Visual comparison: Do observed points fall within simulated bands?</li>
                <li>Calculate goodness-of-fit: Chi-square or Kolmogorov-Smirnov test</li>
                <li>Test for trends: Linear regression on observed vs expected slopes</li>
                <li>Variance analysis: Is observed variance consistent with neutral?</li>
              </ol>
            </div>
            
            <div className="p-4 bg-green-50 rounded">
              <p className="text-sm font-medium mb-2">Key Questions to Answer:</p>
              <ul className="text-sm space-y-1">
                <li>✓ Are operons diverging faster or slower than neutral expectation?</li>
                <li>✓ Do different strains show different patterns?</li>
                <li>✓ Is there evidence for selection (deviation from neutral)?</li>
                <li>✓ Does gene conversion rate match empirical estimates?</li>
              </ul>
            </div>
          </div>
        </div>
      )}

      <div className="bg-gray-100 border-l-4 border-gray-500 p-4 rounded">
        <h3 className="font-semibold mb-2">Study Design Summary:</h3>
        <div className="text-sm space-y-1 text-gray-700">
          <p>• <strong>Null Hypothesis:</strong> Operons evolve neutrally (mutation-drift-gene conversion balance)</p>
          <p>• <strong>Alternative:</strong> Selection acts on rRNA operon similarity/diversity</p>
          <p>• <strong>Test:</strong> Compare observed divergence trajectories to simulated neutral expectations</p>
          <p>• <strong>Parameters from literature:</strong> Gene conversion = 8.6×10⁻⁶, mutation rates from LTEE</p>
        </div>
      </div>
    </div>
  );
};

export default RRNAAnalysis;