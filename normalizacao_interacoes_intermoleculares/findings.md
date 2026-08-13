# Parâmetros de interações intermoleculares

## Resultado

Foi criada uma base inicial para normalização entre o BIOVIA Discovery Studio Visualizer 2024 e o LUNA.

- Discovery Studio: valores recuperados em tempo de execução por `MdmDiscoveryScript`, usando `Mdm::NonbondMonitor::NonbondCriterion()` e `DefaultNonbondCriterion()`.
- LUNA: valores recuperados do `config.cfg` oficial do repositório `keiserlab/LUNA`, versão de código `0.14.0`.
- Nenhum binário ou preferência foi alterado. O único conteúdo novo está nesta pasta.

Os dados completos estão em [dsv_defaults.tsv](./dsv_defaults.tsv) e [luna_defaults.tsv](./luna_defaults.tsv). A primeira correspondência operacional está em [normalized_parameters.tsv](./normalized_parameters.tsv).

## Discovery Studio Visualizer 2024

Instalação consultada: `C:\Program Files\BIOVIA\Discovery Studio 2024\`, versão `24.1.0.23298`.

Principais valores padrão expostos pelo monitor de interações:

| Família | Critério | Valor |
|---|---|---:|
| Ligação de hidrogênio forte | D(H)...A máximo | 3,4 Å |
| Ligação de hidrogênio fraca | D(H)...A máximo | 3,8 Å |
| Ponte salina | D(H)...A máximo | 4,0 Å |
| Ângulos de H-bond | D-H-A, H-A-Y, X-D-A e D-A-Y mínimos | 90° |
| Carga-carga | distância máxima | 5,6 Å |
| Cátion-π / ânion-π | distância máxima | 5,0 Å |
| Cátion-π / ânion-π | ângulo máximo | 40° |
| π-π | centroide máximo | 6,0 Å |
| π-π | átomo mais próximo máximo | 4,5 Å |
| π-π empilhada | theta / gamma máximos | 50° / 35° |
| π-π T-shaped | theta máximo / gamma mínimo | 30° / 55° |
| Halogênio-F | distância máxima | 3,7 Å |
| Halogênio Cl/Br/I | fração VDW máxima | 1,0 |
| Halogênio | C-X-B mínimo / X-B-Y mínimo | 120° / 75° |
| Sulfur-π face-on | distância / ângulo máximos | 4,5 Å / 25° |
| Par solitário-π | distância / ângulo máximos | 3,0 Å / 45° |
| Alkyl | distância de centroide máxima | 5,5 Å |
| VDW bump | fração VDW mínima | 0,7 |

Tipos de interação listados pela API: `saltBridgeType`, `attractiveChargeType`, `waterBridgeHBondType`, `waterHBondType`, `conventionalHBondType`, `carbonHBondType`, `metalAcceptorType`, tipos de halogênio, enxofre, interações π, `alkylType`, `piAlkylType` e tipos desfavoráveis.

## LUNA

O LUNA não está instalado no ambiente Python local consultado. A fonte oficial define os padrões no arquivo `luna/interaction/config.cfg`; a versão atual no código-fonte é `0.14.0`.

Principais valores padrão:

| Família | Parâmetros principais |
|---|---|
| H-bond | D-A ≤ 3,9 Å; H-A ≤ 2,8 Å; D-H-A, H-A-R e D-A-R ≥ 90° |
| H-bond fraca | D-A ≤ 4,0 Å; H-A ≤ 3,0 Å; D-H-A ≥ 110°; variante aromática D-C ≤ 4,5 Å, H-C ≤ 3,5 Å, D-H-C ≥ 120° |
| Iônica | atração e repulsão ≤ 6,0 Å |
| Empilhamento aromático | centroide ≤ 6,0 Å; inclinação 30–60°; deslocamento 30–60° |
| Amida-π | centroide ≤ 4,5 Å; diedro ≤ 30°; deslocamento ≤ 30° |
| Hidrofóbica | ≤ 4,5 Å; superfície mínima de 1 átomo e 1 átomo de interseção |
| Cátion-π | ≤ 6,0 Å |
| Ligação de halogênio | X-A ≤ 4,0 Å; X-centroide ≤ 4,5 Å; C-X-A ≥ 120°; X-A-R ≥ 80°; deslocamento ≤ 60° |
| Ligação de calcogênio | Y-A ≤ 4,0 Å; Y-centroide ≤ 4,5 Å; R-Y-A ≥ 120°; Y-A-N ≥ 80°; deslocamento ≤ 60° |
| Multipolar | N-E ≤ 4,0 Å; ângulos 70–110°; deslocamento ≤ 40° |
| Íon-multipolo | I-D ≤ 4,5 Å; I-D-Y ≥ 60°; deslocamento ≤ 40° |
| Proximal | 2,0–6,0 Å |
| VDW / clash | tolerâncias 0,1 / 0,6 |
| Metal | metal-aceptor ≤ 2,8 Å |

Além dos critérios geométricos, o LUNA possui opções que alteram o conjunto observado: interações covalentes e não covalentes, interações dependentes, pontes de água, regras estritas de doador, separação mínima por ligações e cutoff de região de ligação. Esses campos precisam ser guardados separadamente na normalização; eles não são “parâmetros geométricos” equivalentes aos critérios do DSV.

## Diferenças que impedem uma normalização numérica ingênua

1. **H-bond:** o DSV expõe principalmente D(H)...A e ângulos; o LUNA separa D-A e H-A. Uma linha única perderia informação.
2. **Ponte salina versus interação iônica:** 4,0 Å no DSV é um critério de ponte salina dependente de H-bond; 6,0 Å no LUNA é um limite mais amplo para atração iônica.
3. **π-π:** os dois usam centroide de 6,0 Å, mas classificam orientação por convenções angulares diferentes.
4. **Halogênio e calcogênio/enxofre:** o DSV mistura critérios por elemento e fração VDW; o LUNA usa modelos X-A/Y-A com distâncias absolutas.
5. **VDW/clash:** `0,7` no DSV e `0,6` no LUNA não devem ser tratados como a mesma grandeza sem reimplementar e comparar as fórmulas.
6. **Ausência de parâmetro:** quando uma API não expõe um critério — por exemplo, distância de metal no conjunto consultado do DSV — o valor deve ser armazenado como ausente, não inferido.

## Esquema recomendado para a normalização

Cada registro deve conter:

`software`, `software_version`, `source_kind`, `family`, `native_parameter`, `value`, `unit`, `operator`, `atom_roles`, `interaction_type`, `prerequisites`, `default_or_custom`, `notes`.

Recomendação prática: manter perfis nativos separados (`dsv_2024_default` e `luna_0.14_default`) e criar um terceiro perfil normalizado apenas para campos realmente comparáveis. Isso preserva a auditabilidade e evita afirmar equivalência entre classificadores diferentes.

## Reprodução do levantamento do DSV

O script [probe_dsv_defaults.pl](./probe_dsv_defaults.pl) cria um documento em memória e consulta os defaults expostos pelo monitor. Em PowerShell:

```powershell
$perl = 'C:\Program Files\BIOVIA\Discovery Studio 2024\bin\perl.exe'
$inc = 'C:\Program Files\BIOVIA\Discovery Studio 2024\lib\vendor_perl\5.36.0'
$bin = 'C:\Program Files\BIOVIA\Discovery Studio 2024\bin'
$env:Path = "$bin;$env:Path"
& $perl "-I$inc" "-I$bin" .\probe_dsv_defaults.pl
```

Documentação local relevante:

- [Interaction Criteria dialog](C:/Program%20Files/BIOVIA/Discovery%20Studio%202024/share/doc/DSV/content/docds/client/dl_edit_nonbondmonitor_parameters.htm)
- [AllNonbondCriterionTypes](C:/Program%20Files/BIOVIA/Discovery%20Studio%202024/share/doc/DSV/content/docds/scripting/allnonbondcriteriontypes_function.htm)
- [NonbondCriterion property](C:/Program%20Files/BIOVIA/Discovery%20Studio%202024/share/doc/DSV/content/docds/scripting/nonbondcriterion_property.htm)

Referências públicas utilizadas para o LUNA: [config.cfg oficial](https://raw.githubusercontent.com/keiserlab/LUNA/master/luna/interaction/config.cfg), [documentação do InteractionCalculator](https://luna-toolkit.readthedocs.io/en/latest/api/luna.interaction.calc.html), [código da versão](https://raw.githubusercontent.com/keiserlab/LUNA/master/luna/version.py).
